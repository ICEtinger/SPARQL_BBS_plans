#include "CJTgraph.h"

// CJoin-Table Bipartite Graph
CJTgraph::CJTgraph() {}
	
// initialization  from collection of decomposed triples (i.e.: without ?, *, +, ...)
CJTgraph::CJTgraph(std::vector<triple> q) {
	// m maps variable names (string) to nodes representing them (CJoin*).
	std::unordered_map<std::string, CJoin*> m;
	// for each WHERE triple in query q.
	for (const triple& t : q) {
		Table *table = new Table{ t[1] };
		auto it_kvpair_subject = m.find(t[0]);
		auto it_kvpair_object  = m.find(t[2]);
		CJoin *s = (it_kvpair_subject == m.end() ? new CJoin{ t[0] } : it_kvpair_subject->second);
		CJoin *o = (it_kvpair_object  == m.end() ? new CJoin{ t[2] } : it_kvpair_object->second);
		if (it_kvpair_subject == m.end())
			m[t[0]] = s;
		if (it_kvpair_object == m.end())
			m[t[2]] = o;
		add_edge(table, s, 's'); // 's' is the subject column.
		add_edge(table, o, 'o'); // 'o' is the object column.
	}
}

// initialization from another CJTgraph. 
// Optional parameters m_cj and m_t are maps from the original CJoin/tables to the new ones
// They will be filled by the function if their references are passed as parameters. 
// Else new maps will be created and destroyed within the function. 
CJTgraph::CJTgraph(const CJTgraph& other, std::unordered_map<CJoin*,CJoin*>& m_cj = std::unordered_map<CJoin*, CJoin*>(), std::unordered_map<Table*,Table*>& m_t = std::unordered_map<Table*, Table*>()) {
	// copies g using maps m_cj and m_t.
	for (CJoin* original_cj : other.cjoins) {
		CJoin* cj = new CJoin{ original_cj->name };
		cjoins.insert(cj);
		m_cj.emplace(original_cj, cj);
	}
	for (Table* original_table : other.tables) {
		Table* table = new Table{ original_table->name };
		tables.insert(table);
		m_t.emplace(original_table, table);
	}
	for (CJoin* original_cj : other.cjoins)
		for (const std::pair<Table*, char>& original_edge : original_cj->edges)
			add_edge_without_checking(m_t[original_edge.first], m_cj[original_cj], original_edge.second);
}

CJTgraph::~CJTgraph() {
	for (Table* t : tables)
		delete t;
	for (CJoin* cj : cjoins)
		delete cj;
}

std::unordered_set<Table*>* CJTgraph::get_tables() { return &tables; }
std::unordered_set<CJoin*>* CJTgraph::get_cjoins() { return &cjoins; }

void CJTgraph::print_itself(FILE * pFile = stdout, char format = 'O') {
	if (format == 'O') {
		// prints the operations.
		fprintf(pFile, "n-ary join on:  ");
		std::unordered_set<std::string> joined_tables;
		for (CJoin* cj : cjoins) {
			for (const std::pair<Table*, int> & e : cj->edges)
				joined_tables.insert(e.first->name);
			if (cj->edges.size() >= 2) {
				for (const std::pair<Table*, int> & e : cj->edges)
					fprintf(pFile, "%s.%c=", e.first->name.c_str(), e.second);
				fprintf(pFile, "\b \b ^ "); // "\b \b" deletes last character (in the case, deletes "=").
			}
		}
		fprintf(pFile, "\b \b\b \b"); // "\b \b" deletes last character (in the case, deletes "^ ").
		fprintf(pFile, " on tables: ");
		for (const std::string& table : joined_tables)
			fprintf(pFile, "%s, ", table.c_str());
		fprintf(pFile, "\b \b\b \b\n"); // "\b \b" deletes last character (in the case, deletes ", ").
	}
	else if (format == 'S') {
		// prints the operations.
		fprintf(pFile, "n-ary join on:  ");
		std::unordered_set<std::string> joined_tables;
		for (CJoin* cj : cjoins) {
			for (const std::pair<Table*, int> & e : cj->edges)
				joined_tables.insert(e.first->name);
			if (cj->name[0] == '?')
				fprintf(pFile, "[%s]", cj->name.c_str());
			for (const std::pair<Table*, int> & e : cj->edges)
				fprintf(pFile, "S(%s).%c=", e.first->name.c_str(), e.second);
			fprintf(pFile, "\b \b ^ "); // "\b \b" deletes last character (in the case, deletes "=").
		}
		fprintf(pFile, "\b \b\b \b"); // "\b \b" deletes last character (in the case, deletes "^ ").
		fprintf(pFile, " on tables: ");
		for (const std::string& table : joined_tables)
			fprintf(pFile, "%s, ", table.c_str());
		fprintf(pFile, "\b \b\b \b\n"); // "\b \b" deletes last character (in the case, deletes ", ").
	}
	else {
		// prints in the Graphviz standard format (graphviz.org/)
		// e.g.: digraph G{ "CJoin_0" -> "T_A" [label = a]; "CJoin_0" -> "T_B" [label = b]; "CJoin_0" -> "T_C" [label = c]; }
		// represents a 3-join on the a^{th} column of T_A table with the b^{th} column of T_B with the c^{th} column of T_C.
		fprintf(pFile, "digraph G{ ");
		for (CJoin* cj : cjoins)
			for (const std::pair<Table*, int> & e : cj->edges)
				fprintf(pFile, "\"%s\" -> \"%s\" [label = \"%c\"]; ", cj->name.c_str(), e.first->name.c_str(), e.second);
		fprintf(pFile, " }\n");
	}
}

void CJTgraph::add_edge(Table* t, CJoin* cj, char column) {
	t->edges.insert(cj);
	cj->edges.emplace_back(t, column);
	tables.insert(t);
	cjoins.insert(cj);
}

// does not check for t and cj being on tables and cjoins unordered_sets
void CJTgraph::add_edge_without_checking(Table* t, CJoin* cj, char column) {
	t->edges.insert(cj);
	cj->edges.emplace_back(t, column);
}

// substitute an edge of cj to a new one.
void CJTgraph::substitute_edge(Table* t, CJoin* cj, char column, size_t index_of_edge_to_remove) {
	t->edges.insert(cj);
	cj->edges[index_of_edge_to_remove] = std::make_pair(t, column);
	tables.insert(t);
}

void CJTgraph::evaluate(CJoin* cj) {
	// returns if cjoin outside of the graph.
	if (cjoins.find(cj) == cjoins.end()) {
		printf("ERROR! CJoin outside of this graph");
		return;
	}
	// constructs a vector of the physical tables to join.
	std::vector<Physical_table*> phys_tables_to_join;
	for (const std::pair<Table*, char>& edge : cj->edges)
		phys_tables_to_join.push_back(edge.first->phys_table);

	// constructs a vector of the column positions to be joined in each table. 
	std::vector<char> columns_to_join;
	for (const std::pair<Table*, char>& edge : cj->edges)
		columns_to_join.push_back(edge.second);

	// constructs the name to represent the joined table.
	std::string cjoin_name = cj->name + "(";
	for (const std::pair<Table*, char>& edge : cj->edges)
		cjoin_name += (edge.first->name + "(" + std::to_string(edge.second) + "),");
	cjoin_name[cjoin_name.length() - 1] = ')'; // removes last ','

	// calls function to physically evaluate the join operation.
	std::pair< Physical_table*, std::vector<std::vector<char>> > eval_result = physical_cjoin_evaluator(phys_tables_to_join, columns_to_join);

	Table joined_table{ cjoin_name, std::unordered_set<CJoin*>(), eval_result.first };

	// transfer edges from just-joined tables to the new table.
	// for each table connected to cj:
	for (size_t i = 0, size = cj->edges.size(); i < size; ++i) {
		Table* table = cj->edges[i].first;
		// for each CJoin connected to table:
		for (CJoin* new_neighboring_join : table->edges) {
			// if that CJoin is not the one being evaluated (i.e.: cj)
			if (new_neighboring_join != cj) {
				// then transfer the edge between that CJoin and the just-joined Table to the new Table.
				size_t k = 0;
				// sets k to the index of the edge (between new_neighboring_join and the just-joined table) to be transfered.
				while (new_neighboring_join->edges[k].first != table)
					++k;
				// remove old edge (in index k) and add new one.
				// eval_result.second = M such that M[i][j] is the new column (in the resulted table) of column j of v[i].
				substitute_edge(&joined_table, new_neighboring_join, eval_result.second[i][columns_to_join[i]], k);
			}
		}
		// remove table from graph.
		tables.erase(table);
	}
	// remove cjoin from graph
	cjoins.erase(cj);
}

// explores all original graph evaluation subplans (OGES) and applies function_to_apply to each of them,
// with the variables (OGES, original_graph, cj_using_summaries, edges_indices_ofCJoin, edges_indices_of_Table, edges)
// those 6 variables are created inside the function. First 3 vary between different each OGES. Last 3 are constant.
// parameter q is a query formed only of decomposed triples (i.e.: without ?, *, +, ...)
void CJTgraph::for_each_OGES(const std::vector<triple>& q, std::function< void(CJTgraph&, CJTgraph&, std::unordered_set<CJoin*>&, 
			std::unordered_map<CJoin*,std::vector<size_t>>&, std::unordered_map<Table*,std::vector<size_t>>&, 
			const std::vector<std::tuple<CJoin*,Table*,char>>&, const std::vector<bitArray>&) > function_to_apply) {

	// at first constructs the original graph evaluation subplan that does not use summaries.
	// OGES will be the graph of original-graph evaluation subplans.
	CJTgraph OGES(q);
	std::unordered_map<CJoin*, CJoin*> m_cj;
	std::unordered_map<Table*, Table*> m_t;
	// constructs copy of graph, and m_cj and m_t mapping CJoin and Tables from OGES to original_graph.
	// original_graph will be the graph of summary evaluation subplans.
	CJTgraph original_graph(OGES, m_cj, m_t);
	
	// edges are stored in a vector. There are unordered_maps mapping CJoins and Tables to indices of their respective edges.
	std::unordered_map<CJoin*, std::vector<size_t>> edges_indices_of_CJoin;
	std::unordered_map<Table*, std::vector<size_t>> edges_indices_of_Table;
	std::vector< std::tuple<CJoin*, Table*, char> > edges;
	for (CJoin* cj : *original_graph.get_cjoins()) {
		for (const std::pair<Table*, char>& e : cj->edges) {
			edges_indices_of_CJoin[cj].push_back(edges.size());
			edges_indices_of_Table[e.first].push_back(edges.size());
			edges.emplace_back(cj, e.first, e.second);
		}
	}
	// edges_connections[i] represents all connections of edge i to other edges (by a CJoin or a Table).
	// edges_connections[i] & (1 << j) = (edges i and j are connected).
	std::vector<bitArray> edges_connections(edges.size());
	for (size_t i = 0, size = edges.size(); i < size; ++i) {
		for (const size_t connected_edge_index : edges_indices_of_CJoin[std::get<0>(edges[i])])
			edges_connections[i] |= ((bitArray)1 << connected_edge_index);
		for (const size_t connected_edge_index : edges_indices_of_Table[std::get<1>(edges[i])])
			edges_connections[i] |= ((bitArray)1 << connected_edge_index);
	}

	std::unordered_set<CJoin*> cj_using_summaries;

	// creates nested functions to be called on DFS.

	size_t total_CJoins_added = 0;
	// transforms the graph so that a summary is used on cj.
	auto choose_to_use_summary = [&](CJoin* cj) {
		cj_using_summaries.insert(m_cj[cj]);
		Table* extent_table = new Table{ "E(" + cj->name + ")" };
		// added_cjoins has the cjoins added by this function. They will later be passed to choose_not_to_use_summary to be removed.
		std::vector<CJoin*> added_cjoins;
		added_cjoins.reserve(cj->edges.size());
		// for each table connected to cj, reconnect the table to a newly created cjoin, which is connected to extent_table.
		for (const std::pair<Table*, char>& edge : cj->edges) {
			CJoin* added_cj = new CJoin{ "CJoin_" + std::to_string(total_CJoins_added++) };
			added_cjoins.push_back(added_cj);
			OGES.add_edge_without_checking(extent_table, added_cj, 'e');
			OGES.add_edge_without_checking(edge.first, added_cj, edge.second);
			edge.first->edges.erase(cj);
		}
		OGES.get_tables()->insert(extent_table);
		OGES.get_cjoins()->insert(added_cjoins.begin(), added_cjoins.end());
		OGES.get_cjoins()->erase(cj);
		// cj (and the edges it stores) is still present and unchanged, but inacessible to any table or OGES.cjoins.
		return std::make_pair(extent_table, added_cjoins);
	};

	// "de-transforms" the graph to its form before choose_to_use_summary was applied.
	auto choose_not_to_use_summary = [&](CJoin* cj, std::pair< Table*, std::vector<CJoin*> >& added_t_cj) {
		cj_using_summaries.erase(m_cj[cj]);
		// cj (and the edges it stores) is still present and unchanged after choose_to_use_summary, 
		// but inacessible to any table or OGES.cjoins.
		// now we simply have to add it to OGES.cjoins and make tables point to him as edges again.
		OGES.get_cjoins()->insert(cj);
		for (const std::pair<Table*, char>& edge : cj->edges)
			edge.first->edges.insert(cj);
		// delete and erase added_t and added_cjs.
		total_CJoins_added -= added_t_cj.second.size();
		delete (added_t_cj.first);
		for (CJoin* added_cj : added_t_cj.second)
			delete added_cj;
		OGES.get_tables()->erase(added_t_cj.first);
		for (CJoin* added_cj : added_t_cj.second)
			OGES.get_cjoins()->erase(added_cj);
	};

	// vector of cjoins to be applied DFS on (at first at vector index 0, then 1, ...) 
	std::vector<CJoin*> initial_cjoins(OGES.get_cjoins()->begin(), OGES.get_cjoins()->end());
	// sorts initial_cjoins on the number of edges they have.
	std::sort(initial_cjoins.begin(), initial_cjoins.end(), [](auto cj1, auto cj2){return cj1->edges.size() < cj2->edges.size();});


	// nested function to perform an DFS, and apply function_to_apply at each state of OGES during the exploration.
	std::function<void(size_t, bool)> DFS = [&](size_t node_to_explore, bool apply) {
		if (apply)
			function_to_apply(OGES, original_graph, cj_using_summaries, edges_indices_of_CJoin, edges_indices_of_Table, edges, edges_connections);
		if (node_to_explore == initial_cjoins.size())
			return;
		std::pair< Table*, std::vector<CJoin*> > changed_tables_cjoins = choose_to_use_summary(initial_cjoins[node_to_explore]);
		DFS(node_to_explore + 1, true);
		choose_not_to_use_summary(initial_cjoins[node_to_explore], changed_tables_cjoins);
		DFS(node_to_explore + 1, false);
	};
	
	// now apply the DFS
	DFS(0, true);

}

// explores all original graph evaluation subplans (OGES) and applies function_to_apply for each of them.
// parameter q is a query formed only of decomposed triples (i.e.: without ?, *, +, ...)
std::unordered_set<bitArray> CJTgraph::SES_set(CJTgraph& OGES, CJTgraph& original_graph, std::unordered_set<CJoin*>&  cj_using_summaries,
			std::unordered_map<CJoin*,std::vector<size_t>>& edges_indices_of_CJoin, 
			std::unordered_map<Table*,std::vector<size_t>>& edges_indices_of_Table,
			const std::vector<std::tuple<CJoin*,Table*,char>>& edges, const std::vector<bitArray>& edges_connections) {
	
	// we represent a SES by 64 bits, stored together in a 64-bit int typedefined as "bitArray".
	// to change the number of bits used, it suffices to change the typedef of "bitArray".
	// each bit represents whether an edge is present in the represented SES.
	// the i^th rightest bit (0 <= i < 64) represents the i^th edge on "edges".
	// i.e.: for a bitArray curr_SES, we have curr_SES ^ (1 << i) == true iff edges[i] is present in the represented SES. 

	std::unordered_set<bitArray> SES_visited;

	std::vector<std::pair<bitArray, bitArray>> initial_SES{ {0,0} };
	// adds to initial_SES all SES that have exactly 1 edge connected to each cj having its summary used.
	// for each cj having its summary used:
	for (CJoin* cj : cj_using_summaries) {
		std::vector<std::pair<bitArray, bitArray>> tmp;
		// for each edge connected to that cj:
		for (const size_t edge_index : edges_indices_of_CJoin[cj]) {
			const bitArray edge = (bitArray)1 << edge_index;
			for (const std::pair<bitArray,bitArray>& ses_connections_pair : initial_SES)
				tmp.emplace_back(ses_connections_pair.first | edge, ses_connections_pair.second | edges_connections[edge_index]);
		}
		initial_SES = std::move(tmp);
	}
	for (const std::pair<bitArray, bitArray>& ses_connections_pair : initial_SES)
		SES_visited.insert(ses_connections_pair.first);

	// nested function to be called for Depth First Search.
	std::function<void(bitArray, bitArray)> DFS = [&](bitArray curr_SES, bitArray connected_edges) {
		// all edges that are not present in curr_SES and are connected to curr_SES.
		bitArray possible_next_edges = ~curr_SES & connected_edges;
		// for each edge
		for (size_t i = 0, size = edges.size(); i < size; ++i) {
			const bitArray edge = (bitArray)1 << i;
			// if edge is in possible_next_edges
			if (possible_next_edges & edge) {
				bitArray next_SES = curr_SES | edge;
				// if next_SES has not already been visited, it is inserted in SES_visited and explored by DFS.
				if (SES_visited.insert(next_SES).second)
					DFS(next_SES, connected_edges | edges_connections[i]);
			}
		}
	};

	// run DFS.
	for (const std::pair<bitArray, bitArray>& ses_connections_pair : initial_SES)
		DFS(ses_connections_pair.first, ses_connections_pair.second);

	return SES_visited;
}


int CJTgraph::count_all_OGES(const std::vector<triple>& q) {
	int count = 0;
	for_each_OGES(q, [&count](CJTgraph& OGES, CJTgraph&, std::unordered_set<CJoin*>&,
			std::unordered_map<CJoin*, std::vector<size_t>>&, std::unordered_map<Table*, std::vector<size_t>>&,
			const std::vector<std::tuple<CJoin*, Table*, char>>&, const std::vector<bitArray>&) {++count; });
	return count;
}

void CJTgraph::print_all_OGES(const std::vector<triple>& q, FILE * pFile = stdout, char format = 'O') {
	for_each_OGES(q, [&pFile, &format](CJTgraph& OGES, CJTgraph&, std::unordered_set<CJoin*>&,
			std::unordered_map<CJoin*,std::vector<size_t>>&, std::unordered_map<Table*,std::vector<size_t>>&,
			const std::vector<std::tuple<CJoin*,Table*,char>>&, const std::vector<bitArray>&) {OGES.print_itself(pFile, format); });
}

// explores all original graph evaluation subplans (OGES) and applies function_to_apply for each of them.
// parameter q is a query formed only of decomposed triples (i.e.: without ?, *, +, ...)
void CJTgraph::print_all_SES_for_given_OGES(CJTgraph& OGES, CJTgraph& original_graph, std::unordered_set<CJoin*>&  cj_using_summaries,
			std::unordered_map<CJoin*,std::vector<size_t>>& edges_indices_of_CJoin, 
			std::unordered_map<Table*,std::vector<size_t>>& edges_indices_of_Table, const std::vector<std::tuple<CJoin*,Table*,char>>& edges, 
			const std::vector<bitArray>& edges_connections, FILE * pFile = stdout, char format = 'O') {

	std::unordered_set<bitArray> SES_set_for_given_OGES = SES_set(OGES, original_graph, cj_using_summaries,
			edges_indices_of_CJoin, edges_indices_of_Table, edges, edges_connections);
	for (const bitArray SES : SES_set_for_given_OGES) {
		// clears original graph.
		for (CJoin* cj : *original_graph.get_cjoins())
			cj->edges.clear();
		for (Table* t : *original_graph.get_tables())
			t->edges.clear();
		original_graph.get_cjoins()->clear();
		original_graph.get_tables()->clear();
		// reconstruct only the edges of SES.
		for (size_t i = 0, size = edges.size(); i < size; ++i)
			if (SES & (1 << i))
				original_graph.add_edge(std::get<1>(edges[i]), std::get<0>(edges[i]), std::get<2>(edges[i]));
		original_graph.print_itself(pFile,format);
	}
}

void CJTgraph::print_all_SES_and_OGES(const std::vector<triple>& q, FILE * pFile = stdout, char format = 'O') {
	for_each_OGES(q, [&pFile, &format](CJTgraph& OGES, CJTgraph&original_graph,
			std::unordered_set<CJoin*>& cj_using_summaries, std::unordered_map<CJoin*,std::vector<size_t>>& edges_indices_of_CJoin,
			std::unordered_map<Table*, std::vector<size_t>>& edges_indices_of_Table,
			const std::vector<std::tuple<CJoin*, Table*, char>>& edges, const std::vector<bitArray>& edges_connections) {
		OGES.print_itself(pFile, format);
		fprintf(pFile, "SES:\n");
		print_all_SES_for_given_OGES(OGES, original_graph, cj_using_summaries, edges_indices_of_CJoin, edges_indices_of_Table, 
				edges, edges_connections, pFile, 'S');
		fprintf(pFile, "\n");
	});
}

// returns a vector of decomposed queries. A decomposed query is a vector of triples.
// q is not decomposed. e.g.: a predicate of q is "(A/B)?/C+".
// We decompose that predicate to {C, ABC, CC, ABCC, CCC, ABCCC, ... },
// and queries: { { ? N0 C ? N1}, { ? N0 A ? extra0, ? extra0 B ? extra1, ? extra1 C ? N1 }, { ? N0 C ? extra0, ? extra0 C ? N1 }, ...}.
std::vector<std::vector<triple>> CJTgraph::decompose_query(const std::vector<triple>& q, const int limit_regular_exp) {
	
	// Nested function to be recursively called. Returns the decomposition of a predicate composed by "tokens".
	std::function<std::vector<std::vector<std::string>>(const std::vector<std::string>&)> decompose_predicate = [&] (const std::vector<std::string>& tokens) {
		std::stack<std::string> tokens_read; // except "/", "|".
		std::string binary_operator = ""; // for a binary operator (i.e.: "/" or "|") read.
		// a stack for predicate decompositions. In the end of the function, it will contain only 1 element, which will be returned.
		std::stack<std::vector<std::vector<std::string>>> predicate_decompositions;
		bool any_parenthesis = false; // if there is any "(" in the stack.
		// treats the tokens.
		for (const std::string& token : tokens) {
			// if an "(" was read and not closed yet.
			if (any_parenthesis) {
				if (token == ")") {
					std::stack<std::string> tokens_between_parenthesis;
					while (tokens_read.top() != "(") {
						tokens_between_parenthesis.push(tokens_read.top());
						tokens_read.pop();
					}
					tokens_read.pop(); // pop "(".
					std::vector<std::string> tokens_between_parenthesis_pass;
					while (!tokens_between_parenthesis.empty()) {
						tokens_between_parenthesis_pass.emplace_back(tokens_between_parenthesis.top());
						tokens_between_parenthesis.pop();
					}
					predicate_decompositions.push(decompose_predicate(tokens_between_parenthesis_pass));
				}
				else {
					if (token == "/" || token == "|")
						binary_operator = token;
					else
						tokens_read.push(token);
				}
			}
			// if there is no open "(" read. 
			else {
				// token in { "?", "*", "+", "|" , "/", "(", ")" }, or is a variable name.
				if      (token == "?") { // P? = {(empty expression), P}.
					// adds an empty vector, meaning an empty expression.
					predicate_decompositions.top().push_back(std::vector<std::string>{});
				}
				else if (token == "*") { // P* = {(empty expression), P, PP, PPP, ..., P^limit_regular_exp}.
					for (size_t i = 0, size = predicate_decompositions.top().size(); i < size; ++i) {
						for (int j = 2; j <= limit_regular_exp; ++j) {
							std::vector<std::string> tmp; // tmp will be = P^j.
							for (int k = 0; k < j; ++k) // repeat j times.
								tmp.insert(tmp.end(), predicate_decompositions.top()[i].begin(), predicate_decompositions.top()[i].end());
							predicate_decompositions.top().push_back(std::move(tmp));
						}
					}
					predicate_decompositions.top().push_back(std::vector<std::string>{}); // adds empty vector, meaning empty expression.
				}
				else if (token == "+") { // P+ = {P, PP, PPP, ..., P^limit_regular_exp}.
					for (size_t i = 0, size = predicate_decompositions.top().size(); i < size; ++i) {
						for (int j = 2; j <= limit_regular_exp; ++j) {
							std::vector<std::string> tmp; // tmp will be = P^j.
							for (int k = 0; k < j; ++k) // repeat j times.
								tmp.insert(tmp.end(), predicate_decompositions.top()[i].begin(), predicate_decompositions.top()[i].end());
							predicate_decompositions.top().push_back(std::move(tmp));
						}
					}
				}
				else if (token == "|") { // A|B = {A, B}.
					binary_operator = token;
				}
				else if (token == "/") { // A/B = {AB}. (concatenation).
					binary_operator = token;
				}
				else if (token == "(") {
					tokens_read.push(token);
				}
				else if (token == ")") {
					printf("ERROR: invalid syntax!");
					return std::vector<std::vector<std::string>>();
				} 
				else { // if token is a variable name;
					predicate_decompositions.push(std::vector<std::vector<std::string>>( { {token} } ));
				}
			}
			
			// if there is a binary operator (i.e.: "/" or "|") and two decompositions to evaluate.
			if (binary_operator != "" && predicate_decompositions.size() >= 2) {
				if (predicate_decompositions.size() > 2) {
					printf("ERROR: invalid syntax!");
					return std::vector<std::vector<std::string>>();
				}
				if (binary_operator == "|") {
					std::vector<std::vector<std::string>> popped(std::move(predicate_decompositions.top()));
					predicate_decompositions.pop();
					predicate_decompositions.top().insert(predicate_decompositions.top().end(), popped.begin(), popped.end());
				}
				else { // binary_operator == "/"
					std::vector<std::vector<std::string>> popped_first(std::move(predicate_decompositions.top()));
					predicate_decompositions.pop();
					std::vector<std::vector<std::string>> popped_second(std::move(predicate_decompositions.top()));
					predicate_decompositions.pop();
					std::vector<std::vector<std::string>> to_push;
					for (const std::vector<std::string>& first_elem : popped_second) {
						for (const std::vector<std::string>& second_elem : popped_first) {
							std::vector<std::string> tmp(first_elem); 
							tmp.insert(tmp.end(), second_elem.begin(), second_elem.end());
							to_push.push_back(std::move(tmp));
						}
					}
					predicate_decompositions.push(std::move(to_push));
				}
				binary_operator = "";
			}
		}
		// return the only element left on the predicate_decompositions stack
		return predicate_decompositions.top();
	};
	
	std::vector<std::vector<triple>> total_decomposed{ std::vector<triple>{} };
	
	// for each triple in q, we decompose it into a vector of triples, and "set-multiply" total_decomposed by that vector.
	// e.g.: for triple is "?N0 (A/B)?/C+ ?N1", we have:
	//       triple_decomposition_predicates = {C, ABC, CC, ABCC, CCC, ABCCC, ...}.
	//       triple_decomposition = {{?N0 C ?N1}, {?N0 A ?extra0, ?extra0 B ?extra1, ?extra1 C ?N1}, {?N0 C ?extra0, ?extra0 C ?N1}, ...}.
	for (const triple& to_decompose : q) {
		
		// tokenizes (parses into tokens) the predicate of to_decompose.
		const std::string& predicate = to_decompose[1];
		std::vector<std::string> tokens_parsed;
		std::unordered_set<char> delimiters{ '?', '*', '+', '|' , '/', '(', ')' };
		size_t begin = 0, len = predicate.length();
		while (begin < len) {
			if (delimiters.find(predicate[begin]) != delimiters.end()) {
				tokens_parsed.emplace_back(predicate.substr(begin, 1));
				++begin;
			}
			else {
				size_t end = begin;
				while (delimiters.find(predicate[end]) == delimiters.end() && end < len)
					++end;
				tokens_parsed.emplace_back(predicate.substr(begin, end - begin));
				begin = end;
			}
		}
		
		// the predicates of the triple to_decompose decomposition.
		// e.g.: triple_decomposition_predicates = {C, ABC, CC, ABCC, CCC, ABCCC, ...}
		std::vector<std::vector<std::string>> triple_decomposition_predicates = decompose_predicate(tokens_parsed);

		// now adds extra variables between the decomposed predicates to get the decomposed triples.
		// e.g.: triple_decomposition = {{?N0 C ?N1}, {?N0 A ?extra0, ?extra0 B ?extra1, ?extra1 C ?N1}, {?N0 C ?extra0, ?extra0 C ?N1}, ...}.
		std::vector<triple> a; ///////////////////////////////
		std::vector<std::vector<triple>> triple_decomposition(triple_decomposition_predicates.size(), a);
		// triple_decomposition is a vector of the possible options whose union is the original triple to_decompose.
		// for each such option:
		for (size_t option = 0, s = triple_decomposition.size(); option < s; ++option) {
			// for each predicate in that option:
			for (size_t j = 0, size = triple_decomposition_predicates[option].size(); j < size; ++j)
				triple_decomposition[option].push_back(triple{ "?extra" + std::to_string(j), std::move(triple_decomposition_predicates[option][j]), "?extra" + std::to_string(j + 1) });
			triple_decomposition[option][0][0] = to_decompose[0];
			triple_decomposition[option][triple_decomposition_predicates[option].size() - 1][2] = to_decompose[2];
		}

		// "set-multiply" total_decomposed by the decomposition of triple to_decompose.
		std::vector<std::vector<triple>> tmp;
		for (const std::vector<triple>& decomposed_triple : triple_decomposition) {
			for (const std::vector<triple>& query : total_decomposed) {
				std::vector<triple> new_query(query);
				new_query.insert(new_query.end(), decomposed_triple.begin(), decomposed_triple.end());
				tmp.push_back(std::move(new_query));
			}
		}
		total_decomposed = std::move(tmp);
	}

	return total_decomposed;
}

std::vector<std::vector<triple>> CJTgraph::decompose_query(const std::string& q_to_parse, const int limit_regular_exp) {
	std::vector<triple> parsed_query;
	size_t begin = 0;
	while (q_to_parse[begin] != '{')
		++begin;
	// now we have q_to_parse[begin] = '{'.
	while (true) {
		// now we have q_to_parse[begin] = '{' or '.'

		size_t first_space = begin + 2;
		while (q_to_parse[first_space] != ' ')
			++first_space;
		size_t second_space = first_space + 2;
		while (q_to_parse[second_space] != ' ')
			++second_space;
		size_t third_space = second_space + 2;
		while (q_to_parse[third_space] != ' ' && q_to_parse[third_space] != '.' && q_to_parse[third_space] != '}')
			++third_space;
		// now we have q_to_parse[first_space] = q_to_parse[second_space] = ' ', and q_to_parse[third_space] = ' ' or '.'.
		
		parsed_query.push_back(triple{ q_to_parse.substr(begin + 1, first_space - (begin + 1)), 
		                               q_to_parse.substr(first_space + 1, second_space - (first_space + 1)),
		                               q_to_parse.substr(second_space + 1, third_space - (second_space + 1)) });
		begin = third_space;
		while (q_to_parse[begin] != '.' && q_to_parse[begin] != '}')
			++begin;
		if (q_to_parse[begin] == '}')
			break;
	}
	return decompose_query(parsed_query, limit_regular_exp);
}

void CJTgraph::print_decomposed_query(const std::vector<triple>& decomposed_query, FILE* pFile) {
	fprintf(pFile, "decomposed query: SELECT (...) WHERE{\n");
	for (const triple& tr : decomposed_query)
		fprintf(pFile, "   %s %s %s\n", tr[0].c_str(), tr[1].c_str(), tr[2].c_str());
	fprintf(pFile, "}\n");
}

int main() {
	printf("STARTED RUNNING! \n");
	std::vector<triple> q{ triple{"?N0", "A","?N1" }, triple{ "?N1", "B","?N2" } };
	//std::vector<triple> q{ triple{ "?N0", "A+","?N1" } };

	FILE * pFile = stdout; 
	//FILE * pFile = fopen("output_filename", "w");
	
	//print_all_OGES(q, pFile);
	CJTgraph::print_all_OGES(q, pFile);

	//CJTgraph::print_all_SES_and_OGES(q, pFile);
	
	std::vector<std::vector<triple>> q_decomposition = CJTgraph::decompose_query(q, 3);

	/*
	for (const std::vector<triple>& decomposed_q : q_decomposition) {
		CJTgraph::print_decomposed_query(decomposed_q, pFile);
		fprintf(pFile, "OGES:\n");
		CJTgraph::print_all_OGES(decomposed_q, pFile);
		fprintf(pFile, "\n\n");
	}
	*/

	// https://en.wikipedia.org/wiki/Shunting-yard_algorithm
	
	int wait;
	scanf("%d", &wait);

	return 0;
}

