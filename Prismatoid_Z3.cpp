#include "Prismatoid_Z3.h"

#define DEBUG 0
 
using namespace ogdf; 
using namespace std;

int strToInt(string a)
{
        istringstream iss(a);
        int tmp;
        iss>>tmp;
        return tmp;
}


template<typename type> string toString(type x)
{
	std::ostringstream os;
	os << x;
	return os.str();
}

Prismatoid_Z3::Prismatoid_Z3(std::string file_name, int dim, bool bit_vector):
faces(2), dim(dim), bit_vector(bit_vector)//, nFacet(2, 0)
{
	Z3_config cfg = Z3_mk_config();
	//Z3_set_param_value(cfg, "MODEL", "true");
	ctx = Z3_mk_context(cfg);
	Z3_del_config(cfg);

	Z3_parse_smtlib_file(ctx, file_name.c_str(), 0, 0, 0, 0, 0, 0);

	int num_formulas = Z3_get_smtlib_num_formulas(ctx);
	for (int i = 0; i < num_formulas; i++) {
		Z3_ast f = Z3_get_smtlib_formula(ctx, i);
		//printf("formula %d: %s\n", i, Z3_ast_to_string(ctx, f));
		Z3_assert_cnstr(ctx, f);
		/*if(Z3_check(ctx) == Z3_L_FALSE)
		cout<<"not check"<<endl;*/
	}

	int num_decls = Z3_get_smtlib_num_decls(ctx);
	for (int i = 0; i < num_decls; i++) 
	{
		Z3_func_decl d = Z3_get_smtlib_decl(ctx, i);

		if(d == 0) continue;
		//cout<<Z3_func_decl_to_string(ctx, d);
		if(string(Z3_func_decl_to_string(ctx, d)).find("sk_hack") == string::npos)
		{
			string s(decl_to_str_var(d));
			faces[get_dim(s)].push_back(d);

			/*if(get_dim(s) == 0)
			{
				vector<string> F(2);
				get_facets(faces[0].size() - 1, F);
				for(int j = 0; j < 2; j++)
				{
					for(int k = 0; k < F[j].size(); k++)
					{	
						string s = "";
						s+=F[j][k];		
						if(strToInt(s) + 1 > nFacet[j])
							nFacet[j] = strToInt(s) + 1;
					}
				}
			}*/
		}
	}

	// Fill adj
	adj.resize(faces[1].size());
	adj_e.resize(faces[0].size());
	for(int i = 0; i < faces[1].size(); i++)
	{
		for(int j = 0; j < faces[0].size(); j++)
		{
			string v = decl_to_str_var(faces[0][j]), e = decl_to_str_var(faces[1][i]);
			if(isAdj(v, e))
			{
				adj[i].push_back(j);
				adj_e[j].push_back(i);
			}
		}
	}
}

int Prismatoid_Z3::get_dim(std::string &var)
{
	return dim + 2 - get_nFacets(var);
}

int Prismatoid_Z3::get_nFacets(std::string &var)
{
	return var.size() - 2;
}

std::string Prismatoid_Z3::decl_to_str_var(Z3_func_decl &d)
{
	Z3_ast a = Z3_mk_app(ctx, d, 0, 0);
	return string(Z3_ast_to_string(ctx, a));
}

int Prismatoid_Z3::belong_to(Z3_func_decl &face)
{
	string f = decl_to_str_var(face);
	vector<string> F(2);
	get_facets(f, F);
	if(F[0].size() == 1)
		return 2;
	if(F[1].size() == 1)
		return 1;
	return 0;
}

ogdf::String Prismatoid_Z3::get_color(int vertex)
{
	switch(belong_to(faces[0][vertex]))
	{
	case 1:
		return "#0000FF";
	case 2:
		return "#FF0000";
	}
	return "#000000";
}

void Prismatoid_Z3::get_facets(int var, std::vector<std::string> &F)
{
	string s = decl_to_str_var(faces[0][var]);
	return get_facets(s, F);
}

void Prismatoid_Z3::get_facets(std::string &var, vector<string> &F)
{
	int pos = var.find("x");
	F[0] = var.substr(1, pos - 1);
	F[1] = var.substr(pos + 1, var.npos);
}

bool Prismatoid_Z3::isSubset(std::string &s1, std::string &s2)
{
	for(int i = 0; i < s1.size(); i++)
		if(s2.find(s1[i]) == s2.npos)
			return false;
	return true;
}

bool Prismatoid_Z3::isAdj(string &v, string &e)
{
	vector<string> f_v(2), f_e(2);
	get_facets(v, f_v);
	get_facets(e, f_e);
	return isSubset(f_e[0], f_v[0]) && isSubset(f_e[1], f_v[1]);
}

int Prismatoid_Z3::get_val(Z3_model &m, Z3_func_decl &d)
{
	int ret = 0;	
	/*if(bit_vector)
	{
		Z3_ast v;
		Z3_eval_func_decl(ctx, m, d, &v);
		Z3_ast v_i = Z3_mk_bv2int(ctx, v, false);	
		Z3_get_numeral_int(ctx, v_i, &ret);
	}
	else
	{*/
	Z3_ast v = 0;
	Z3_eval_func_decl(ctx, m, d, &v);
	assert(v != 0);
	Z3_get_numeral_int(ctx, v, &ret);
	//}
	return ret;
}

string Prismatoid_Z3::inter(std::string &s1, std::string &s2)
{
	string s;
	for(int i = 0; i < s1.size(); i++)
		if(s2.find(s1[i]) != s2.npos)
			s += s1[i];
	return s;
}

void Prismatoid_Z3::make_basis(Graph &G, Z3_model &m, int nBasis) 
{
	//vector<int> adj_edge(faces.size(), -1);

	vector<node> allNodes(faces[0].size(), NULL);
	for(int i = 0; i < faces[0].size(); i++)
	{
		if(get_val(m, faces[0][i]) == 0 || (belong_to(faces[0][i]) != nBasis)) continue;
		allNodes[i] = G.newNode(i);
	}	

	for(int i = 0; i < allNodes.size(); i++)
	{
		if(allNodes[i] == NULL) continue;
		for(int j = i+1; j < allNodes.size(); j++)
		{
			if(allNodes[j] == NULL) continue;

			vector<string> F1(2), F2(2);
			get_facets(i, F1);
			get_facets(j, F2);
			string s = inter(F1[nBasis-1], F2[nBasis-1]);
			if(s.size() == dim) // if the vertices are adjacents
				G.newEdge(allNodes[i], allNodes[j]);
		}
	}		
}

void Prismatoid_Z3::make_graph(Graph &G, Z3_model &m)
{
	vector<node> allNodes(faces[0].size(), 0);
	for(int i = 0; i < faces[0].size(); i++)
	{
		if(get_val(m, faces[0][i]) == 0) continue;
		allNodes[i] = G.newNode(i);
	}

	for(int i = 0; i < faces[1].size(); i++)
	{
		if(get_val(m, faces[1][i]) == 0) continue;

		//cout<<"edge "<<Z3_func_decl_to_string(ctx, faces[1][i]);

		int u = -1, v = -1; // first and second extremities
		int j;
		for(j = 0; j < adj[i].size(); j++)
		{
			if(get_val(m, faces[0][adj[i][j]]) == 1) // if adj[i][j] is a "real" vertex
			{
				u = adj[i][j];
				break;
			}
		}
		assert(u != -1);
		for(++j; j < adj[i].size(); j++)
		{
			if(get_val(m, faces[0][adj[i][j]]) == 1) // if adj[i][j] is a "real" vertex
			{
				assert(v == -1);
				v = adj[i][j];
			}
		}
		assert(v != -1);
		node node_u = allNodes[u];
		node node_v = allNodes[v];
		G.newEdge(node_u, node_v, i);
	}
}

void Prismatoid_Z3::write_graph(Graph &G, string file_name)
{
	GraphAttributes GA(G, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |
		GraphAttributes::nodeLabel | GraphAttributes::nodeColor | 
		GraphAttributes::edgeColor | GraphAttributes::edgeStyle | 
		GraphAttributes::nodeStyle | GraphAttributes::nodeTemplate);
	GA.directed(false);
	node v;
	forall_nodes(v,G)
	{
		GA.width(v) = GA.height(v) = 40.0;
		GA.colorNode(v) = get_color(v->index());
	}

	FMMMLayout fmmm;

	fmmm.useHighLevelOptions(true);
	fmmm.unitEdgeLength(550.0); 
	fmmm.newInitialPlacement(true);
	fmmm.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);
	fmmm.call(GA);
	GA.writeGML(file_name.c_str());
}

void Prismatoid_Z3::write_smt(std::string file_name, Z3_model &m)
{
	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	
	string decl = ":extrafuns (";
	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < faces[i].size(); j++)
		{
			decl += " (" + decl_to_str_var(faces[i][j]) + " BitVec[1])";
		}
	}
	outfile<<decl<<" )\n";

	for(int i = 0; i < 2; i++)
	{
		for(int j = 0; j < faces[i].size(); j++)
		{
			outfile<<":formula (= "<<decl_to_str_var(faces[i][j])<<" bv"<<toString(get_val(m, faces[i][j]))<<"[1])\n";
		}
	}
	outfile.close();
}

bool Prismatoid_Z3::bfs_in_face(vector<node> &v, vector< vector<node> > &cut, ogdf::Graph &G, bool enforce)
{
	vector< vector<string> > F(2);
	vector<string> F_inter(2);
	for(int i = 0; i < 2; i++)
	{
		F[i].resize(2);
		get_facets(v[i]->index(), F[i]);
	}
	for(int i = 0; i < 2; i++)
		F_inter[i] = inter(F[0][i], F[1][i]);
#if DEBUG
	cout<<"("<<F[0][0]<<", "<<F[0][1]<<") to ("<<F[1][0]<<", "<<F[1][1]<<")"<<endl;
	cout<<"inter: ("<<F_inter[0]<<", "<<F_inter[1]<<")"<<endl;
#endif
	cut.resize(2);

	for(int i = 0; i < 2; i++) // BFS from v[i]
	{
		vector<bool> visited(faces[0].size(), false);
		std::queue<node> q;
		q.push(v[i]);
		visited[v[i]->index()] = true;
		while(q.size())
		{
			node next = q.front();
			q.pop();
			cut[i].push_back(next);
			if(next == v[1-i])
			{
				return true;
			}
			edge e_adj;
			forall_adj_edges(e_adj, next)
			{
				node adj = e_adj->opposite(next);
				int index = adj->index();
				vector<string> F_adj(2);
				get_facets(index, F_adj);
#if DEBUG
				cout<<"adj: ("<<F_adj[0]<<", "<<F_adj[1]<<")"<<endl;
#endif
				if(!visited[index] && isSubset(F_inter[0], F_adj[0]) && isSubset(F_inter[1], F_adj[1]))
				{
					q.push(adj);
					visited[index] = true;
				}
#if DEBUG
				else
				{
					if(visited[index])
						cout<<"Visited"<<endl;
					else
						cout<<"Not in inter"<<endl;
				}
#endif
			}
		}
		//cout<<"endpoint not found"<<endl<<endl;
	}

	if(enforce)
		ensure_convexity(v[0]->index(), v[1]->index(), F_inter, cut, G);
	return false;
}

Z3_ast Prismatoid_Z3::toBool(Z3_ast a, bool b)
{
	Z3_ast ret;
	if(!bit_vector)
	{
		Z3_sort int_sort = Z3_mk_int_sort(ctx);
		ret = Z3_mk_eq(ctx, a, Z3_mk_numeral(ctx, b ? "1" : "0", int_sort));
	}
	else
	{
		Z3_sort bv_sort = Z3_mk_bv_sort(ctx, 1);
		ret = Z3_mk_eq(ctx, a, Z3_mk_numeral(ctx, b ? "1" : "0", bv_sort));
	}
	return ret;
}

void Prismatoid_Z3::ensure_convexity(int s, int t, std::vector<std::string> &F_inter, std::vector< std::vector<ogdf::node> > &cut, ogdf::Graph &G)
{
	vector<int> nodesNeighborhood; // nodesNeighborhood will contain all the nodes potentially adjacent to a node in cut[0]
	// but not here
	for(int i = 0; i < cut[0].size(); i++)
	{
		int u = cut[0][i]->index();
		for(int j = 0; j < adj_e[u].size(); j++) // for all edges uv with u as a vertex
		{
			int uv = adj_e[u][j];
			for(int k = 0; k < adj[uv].size(); k++) // for all vertices with this edge as a vertex
			{
				int v = adj[uv][k];
				//cout<<" "<<v<<endl;
				vector<string> F_adj(2);
				get_facets(v, F_adj);
				if(!isSubset(F_inter[0], F_adj[0]) || !isSubset(F_inter[1], F_adj[1])) // if v is in the "good" facets
					continue;

				for(int l = 0; l < cut[0].size(); l++)
					if(v == cut[0][l]->index()) // we don't want to add vertices already in G
						goto cont;

				if(find(nodesNeighborhood.begin(), nodesNeighborhood.end(), v) == nodesNeighborhood.end()) 
					// To not add twice v in nodesNeighborhood
					nodesNeighborhood.push_back(v);

cont: continue;
			}
		}
	}

	// add constraints
	//vector<Z3_ast> not_vert(2);
	Z3_ast out_vert;
	/*Z3_ast *all;*/

	//for(int i = 0; i < 2; i++) // If all nodes in one side are not here, we may have a solution
	//{
	//	all = new Z3_ast[cut[i].size()];
	//	for(int j = 0; j < cut[i].size(); j++)
	//	{
	//		int vertex = cut[i][j]->index();
	//		Z3_ast a = Z3_mk_app(ctx, faces[0][vertex], 0, 0);
	//		all[j] = Z3_mk_not(ctx, toBool(a));
	//	}
	//	not_vert[i] = Z3_mk_and(ctx, cut[i].size(), all);
	//	delete[] all;
	//}

	Z3_ast *all_neighbors = new Z3_ast[nodesNeighborhood.size()];
	for(int i = 0; i < nodesNeighborhood.size(); i++)
	{
		all_neighbors[i] = toBool(Z3_mk_app(ctx, faces[0][nodesNeighborhood[i]], 0, 0));
	}
	out_vert = Z3_mk_or(ctx, nodesNeighborhood.size(), all_neighbors);
	delete[] all_neighbors;

	/*Z3_ast all_cond[3];
	for(int i = 0; i < 2; i++)
		all_cond[i] = not_vert[i];
	all_cond[2] = out_vert;*/

	Z3_ast all_cond[3];
	all_cond[0] = toBool(Z3_mk_app(ctx, faces[0][s], 0, 0), false);
	all_cond[1] = toBool(Z3_mk_app(ctx, faces[0][t], 0, 0), false);
	all_cond[2] = out_vert;

	Z3_ast final = Z3_mk_or(ctx, 3, all_cond);
	Z3_assert_cnstr(ctx, final);
	//string s(Z3_ast_to_string(ctx, final));
	//cout<<"final: "<<s<<endl;
}

bool Prismatoid_Z3::check_enforce_convexity(ogdf::Graph &G, bool enforce)
{
	vector<node> allNodes;

	node u;
	forall_nodes(u, G)
		allNodes.push_back(u);

	vector<node> v(2);
	int nbFound = 0;
	for(int i = 0; i < allNodes.size(); i++)
	{
		v[0] = allNodes[i];
		for(int j = i+1; j < allNodes.size(); j++)
		{
			v[1] = allNodes[j];
			vector< vector<node> > cut(2);
			if(!bfs_in_face(v, cut, G, enforce))
			{
				return false;
			}
		}
	}
	return true;
}

int Prismatoid_Z3::check_enforce_kconnectivity(ogdf::Graph &G)
{
	EdgeArray<double> w(G, 1);
	MinCut minCut(G, w);
	return minCut.minimumCut();
}

void Prismatoid_Z3::run(int min_connectivity, bool find_all, bool enforce)
{
	//cout<<toString(nFacet[0])<<"x"<<toString(nFacet[1])<<endl;
	int nFound = 0;
	int nPossible = 0;
	
	string file_graph = "G";// + toString(nFacet[0]) + "x" + toString(nFacet[1]) + "n";
	string file_smt = "S";// + toString(nFacet[0]) + "x" + toString(nFacet[1]) + "n";
	time_t t_d = clock();
	while(true)
	{
		Z3_model m      = 0;

		//Z3_open_log(ctx, "z3.log");
		//Z3_trace_to_file(ctx, "z3.trace");

		time_t t1 = clock();
		Z3_lbool result = Z3_check_and_get_model(ctx, &m);
		time_t t2 = clock();
		//cosut<<"Time " + toString(t2 - t1)<<endl;
		nPossible++;
		if(nPossible % 100 == 0)
			cout<<"Possible: "<<nPossible<<endl;
		
		switch (result) 
		{
		case Z3_L_FALSE:
			printf("unsat\n");
			return;
		case Z3_L_UNDEF:
			printf("unknown\n");
			printf("potential model:\n%s\n", Z3_model_to_string(ctx, m));
			return;
		case Z3_L_TRUE:
			//cout<<"sat\n\n, "<<Z3_model_to_string(ctx, m);
			break;
		}
		//cout<<"found"<<endl;
		//cout<<Z3_model_to_string(ctx, m)<<endl;
		Graph G;

		make_graph(G, m);

		bool convexity = check_enforce_convexity(G, enforce);

		if(!convexity)
		{
			if(!enforce)
			{
				Z3_literals lits = Z3_get_relevant_literals(ctx);
				if (lits) {
					Z3_block_literals(ctx, lits);
					Z3_del_literals(ctx, lits);
				}
			}
			
			//cout<<"Found not convex graph"<<endl;
		}
		else
		{
			if(min_connectivity > 0)
			{
				if(check_enforce_kconnectivity(G) < min_connectivity)
				{
					Z3_literals lits = Z3_get_relevant_literals(ctx);
					if (lits) {
						Z3_block_literals(ctx, lits);
						Z3_del_literals(ctx, lits);
					}
					cout<<"Found convex but not connected"<<endl;
					continue;
				}
			}
			cout<<"time: "<<(t2-t_d)/1000.<<endl;
			cout<<"Found "<<min_connectivity<<"-connected and convex graph"<<endl;

			write_graph(G, file_graph + toString(nFound) + ".gml");
			write_smt(file_smt + toString(nFound) + ".smt", m);

			Graph G1, G2;
			make_basis(G1, m, 1);
			write_graph(G1, file_graph + toString(nFound) + "+.gml");
			
			make_basis(G2, m, 2);
			write_graph(G2, file_graph + toString(nFound) + "-.gml");

			if(!find_all)
				return;
			nFound++;
			Z3_literals lits = Z3_get_relevant_literals(ctx);
			if (lits) {
				Z3_block_literals(ctx, lits);
				Z3_del_literals(ctx, lits);
			}
		}
		//vector<int> vertices_true;
		//vector<int> all_to_true(faces[0].size(), -1);
		//for(int i = 0; i < faces[0].size(); i++)
		//{
		//	if(!get_val(m, faces[0][i])) continue;
		//	//cout<<decl_to_str_var(faces[0][i])<<" ";
		//	all_to_true[i] = vertices_true.size();
		//	vertices_true.push_back(i);
		//}

		//for(int i = 0; i < faces[1].size(); i++)
		//{
		//	if(!get_val(m, faces[1][i])) continue;
		//	
		//	int u = -1, v = -1; // first and second extremities
		//	int j;
		//	for(j = 0; j < adj[i].size(); j++)
		//		if(all_to_true[adj[i][j]] != -1) // if adj[i][j] is a "real" vertex
		//		{
		//			//cout<<decl_to_str_var(faces[0][adj[i][j]])<<" (edge "<<decl_to_str_var(faces[1][i])<<") ";
		//			u = all_to_true[adj[i][j]];
		//			break;
		//		}
		//	assert(u != -1);
		//	for(++j; j < adj[i].size(); j++)
		//		if(all_to_true[adj[i][j]] != -1) // if adj[i][j] is a "real" vertex
		//		{
		//			assert(v==-1);
		//			//cout<<decl_to_str_var(faces[0][adj[i][j]])<<" ";
		//			v = all_to_true[adj[i][j]];
		//			// TODO: break;
		//		}
		//	assert(v != -1);
		//	
		//}


	}
	time_t t_e = clock();
	cout<<"Total time " + toString((t_e - t_d)/1000.)<<endl;
}
