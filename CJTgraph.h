#pragma once

#include <cstdio>
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stack>

// each triple is (subject, predicate, object)
typedef std::array<std::string, 3> triple;
// represents an array of 64 bits.
typedef int_fast64_t bitArray;

struct Physical_table {
	// e.g.: std::vector<std::tuple<>> data; int statistics;
};

// first element of returned pair: a pointer to the resulted Physical_table
// second element: vector<vector<char> M such that M[i][j] is the new column position (in the resulted table) of column j of v[i].
std::pair<Physical_table*, std::vector<std::vector<char>> > physical_cjoin_evaluator(std::vector<Physical_table*> v, std::vector<char> positions_columns_to_join) {
	return std::pair<Physical_table*, std::vector<std::vector<char>> >();
	// ...
}

struct CJoin;
struct Table {
	std::string name;
	std::unordered_set<CJoin*> edges;
	Physical_table* phys_table;
};

struct CJoin {
	std::string name;
	std::vector< std::pair<Table*, char> > edges;
};

// CJoin-Table Bipartite Graph
class CJTgraph {

private:
	std::unordered_set<Table*> tables;
	std::unordered_set<CJoin*> cjoins;

public:
	CJTgraph();

	// initialization  from collection of decomposed triples (i.e.: without ?, *, +, ...)
	CJTgraph(std::vector<triple> q);

	// initialization from another CJTgraph. 
	// Optional parameters m_cj and m_t are maps from the original CJoin/tables to the new ones
	// They will be filled by the function if their references are passed as parameters. 
	// Else new maps will be created and destroyed within the function. 
	CJTgraph(const CJTgraph& other, std::unordered_map<CJoin*, CJoin*>& m_cj, std::unordered_map<Table*, Table*>& m_t);

	~CJTgraph();

	std::unordered_set<Table*>* get_tables();
	std::unordered_set<CJoin*>* get_cjoins();

	// output_filename = empty string to use standard output.
	void print_itself(FILE * pFile, char format);

	void add_edge(Table* t, CJoin* cj, char column);

	// does not check for t and cj being on tables and cjoins unordered_sets
	void add_edge_without_checking(Table* t, CJoin* cj, char column);

	// substitute an edge of cj to a new one.
	void substitute_edge(Table* t, CJoin* cj, char column, size_t index_of_edge_to_remove);

	void evaluate(CJoin* cj);

	// explores all original graph evaluation subplans (OGES) and applies function_to_apply to each of them,
	// with the variables (OGES, original_graph, cj_using_summaries, edges_indices_ofCJoin, edges_indices_of_Table, edges)
	// those 6 variables are created inside the function. First 3 vary between different each OGES. Last 3 are constant.
	// parameter q is a query formed only of decomposed triples (i.e.: without ?, *, +, ...)
	static void for_each_OGES(const std::vector<triple>& q, std::function< void(CJTgraph&, CJTgraph&, std::unordered_set<CJoin*>&,
		std::unordered_map<CJoin*, std::vector<size_t>>&, std::unordered_map<Table*, std::vector<size_t>>&,
		const std::vector<std::tuple<CJoin*, Table*, char>>&, const std::vector<bitArray>&) > function_to_apply);

	// explores all original graph evaluation subplans (OGES) and applies function_to_apply for each of them.
	// parameter q is a query formed only of decomposed triples (i.e.: without ?, *, +, ...)
	static std::unordered_set<bitArray> CJTgraph::SES_set(CJTgraph& OGES, CJTgraph& original_graph, std::unordered_set<CJoin*>&  cj_using_summaries,
		std::unordered_map<CJoin*, std::vector<size_t>>& edges_indices_of_CJoin,
		std::unordered_map<Table*, std::vector<size_t>>& edges_indices_of_Table,
		const std::vector<std::tuple<CJoin*, Table*, char>>& edges, const std::vector<bitArray>& edges_connections);

	static int count_all_OGES(const std::vector<triple>& q);

	static void print_all_OGES(const std::vector<triple>& q, FILE * pFile, char format);

	static void print_all_SES_for_given_OGES(CJTgraph& OGES, CJTgraph& original_graph, std::unordered_set<CJoin*>&  cj_using_summaries,
		std::unordered_map<CJoin*, std::vector<size_t>>& edges_indices_of_CJoin,
		std::unordered_map<Table*, std::vector<size_t>>& edges_indices_of_Table, const std::vector<std::tuple<CJoin*, Table*, char>>& edges,
		const std::vector<bitArray>& edges_connections, FILE * pFile, char format);

	static void print_all_SES_and_OGES(const std::vector<triple>& q, FILE * pFile, char format);

	// returns a vector of decomposed queries. A decomposed query is a vector of triples.
	// q is not decomposed. e.g.: a predicate of q is "(A/B)?/C+".
	// We decompose that predicate to {C, ABC, CC, ABCC, CCC, ABCCC, ... },
	// and queries: { { ? N0 C ? N1}, { ? N0 A ? extra0, ? extra0 B ? extra1, ? extra1 C ? N1 }, { ? N0 C ? extra0, ? extra0 C ? N1 }, ...}.
	static std::vector<std::vector<triple>> decompose_query(const std::vector<triple>& q, const int limit_regular_exp);

	static std::vector<std::vector<triple>> decompose_query(const std::string& q_to_parse, const int limit_regular_exp);

	static void print_decomposed_query(const std::vector<triple>& decomposed_query, FILE* pFile);
};

