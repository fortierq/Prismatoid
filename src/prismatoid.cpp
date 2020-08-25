#include "prismatoid.h"

using namespace std;

template<typename type> string toString(type x)
{
	std::ostringstream os;
	os << x;
	return os.str();
}

Prismatoid::Prismatoid(bool bit_vector, int dim, vector<int> nFacet, string file_name, bool verbose):
dim(dim), file_name(file_name), verbose(verbose), nFacet(nFacet), nGeo(nFacet.size()), bit_vector(bit_vector)
{
	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary); 

	maxFacets = *max_element(nFacet.begin(), nFacet.end());
	outfile<<"(benchmark test\n"; 

	outfile.close();
}

string Prismatoid::neg(std::string s)
{
	return "(not " + s + ")";
}

string Prismatoid::toVariable(vector<Set> &s)
{
	string str = "F";
	for(int i = 0; i < s.size(); i++)
	{
		str += setToString(s[i], nFacet[i]);
		if(i != s.size() - 1)
			str += "x";
	}
	return str;
}

string Prismatoid::toBoolVar(vector<Set> &s, bool b)
{
	string str = "(= ";
	str += toVariable(s);
	if(bit_vector)
	{
		if(b)
			str += " bv1[1])";
		else
			str += " bv0[1])";
	}
	else
	{
		if(b)
			str += " 1)";
		else 

			str += " 0)";
	}	
	return str;
}

string Prismatoid::setToString(Set &s, int d)
{
	string c = "";
	for(int i = 0; i < d; i++)
		if(s.contains(i))
			c += toString(i);
	return c;
}

string Prismatoid::formula(std::string s)
{
	return ":formula (" + s + ")\n";
}

void Prismatoid::assert_0_or_1(int d)
{
	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	outfile<<";\tCONDITION 0 or 1 "<<d<<"-FACES\n";

	vector<Set> s(nGeo);
	for(int k = 1; k < (dim + 2 - d); k++) // dim + 2 - d is the number of facets containing a d-face
		for(s[0] = Set::first_n(k); !s[0].contains(nFacet[0]); s[0] = s[0].colex_next())
			for(s[1] = Set::first_n(dim + 2 - d - k); !s[1].contains(nFacet[1]); s[1] = s[1].colex_next())
				if(!bit_vector)
					outfile<<formula("or (= " + toVariable(s) + " 0) (= " + toVariable(s) + " 1)");
				else
					outfile<<formula("or (= " + toVariable(s) + " bv0[1]) (= " + toVariable(s) + " bv1[1])");
}

void Prismatoid::declare_face(int d)
{
	time_t t1 = clock();

	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	outfile<<";\tDECLARE "<<d<<"-FACES\n";

	outfile<<":extrafuns (";
	vector<Set> s(nGeo);
	for(int k = 1; k < (dim + 2 - d); k++) // dim + 2 - d is the number of facets containing a d-face
		for(s[0] = Set::first_n(k); !s[0].contains(nFacet[0]); s[0] = s[0].colex_next())
			for(s[1] = Set::first_n(dim + 2 - d - k); !s[1].contains(nFacet[1]); s[1] = s[1].colex_next())
				outfile<<" (" + toVariable(s) + (bit_vector ? " BitVec[1])" : " Int)");
	outfile<<" )\n";
	time_t t2 = clock();
	outfile.close();
	if(verbose)
		cout<<d<<"-faces declared"<<endl<<"Time: " + toString(t2 - t1)<<endl;
}

void Prismatoid::assert_no_dstep()
{
	time_t t1 = clock();

	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	outfile<<";\tSEPARATION PROPERTY\n";

	vector< vector<Set> > s(nGeo); // "complete" (with d+1 facets) and other component of the ith vertex 
	s[0].resize(2);
	s[1].resize(2);

	string cond = "(and";
	for(s[0][0] = Set::first_n(dim + 1); !s[0][0].contains(nFacet[0]); s[0][0] = s[0][0].colex_next()) // for all blue vertices
		for(s[1][1] = Set::first_n(dim + 1); !s[1][1].contains(nFacet[1]); s[1][1] = s[1][1].colex_next()) // for all red vertices
			for(s[0][1] = Set::first_n(1); !s[0][1].contains(nFacet[1]); s[0][1] = s[0][1].colex_next())
			{
				if(!s[0][1].subset(s[1][1]))
					continue;
				for(s[1][0] = Set::first_n(1); !s[1][0].contains(nFacet[0]); s[1][0] = s[1][0].colex_next())
				{
					if(!s[1][0].subset(s[0][0]))
						continue;

					//outfile<<formula(neg("(and " + toBoolVar(s[0]) + " " + toBoolVar(s[1]) + ")"));
					cond += " ";
					cond += neg("(and " + toBoolVar(s[0]) + " " + toBoolVar(s[1]) + ")");
				}
			}
			cond += ")";
			outfile<<formula(cond);
			time_t t2 = clock();
			if(verbose)
				cout<<"No d-step property assumed"<<endl<<"Time: " + toString(t2 - t1)<<endl;
}

// warning: minkowski non simple!!!!
void Prismatoid::assert_degree()
{
	time_t t1 = clock();

	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	outfile<<";\tDEGREE CONDITION\n";

	vector<Set> s(nGeo);
	for(int k = 1; k < dim + 1; k++)
		for(s[0] = Set::first_n(k); !s[0].contains(nFacet[0]); s[0] = s[0].colex_next())
			for(s[1] = Set::first_n(dim + 2 - k); !s[1].contains(nFacet[1]); s[1] = s[1].colex_next())
			{
				string adjEdges = "(and";
				for(int j = 0; j < 2; j++)
				{
					if(s[j].cardinality() == 1) continue;
					for(int i = 0; i < nFacet[j]; i++)
					{
						if(!s[j].contains(i)) continue; 
						vector<Set> tmp = s;
						tmp[j] = s[j] - i;
						adjEdges += (" " + toBoolVar(tmp));
					}
				}
				adjEdges += ")";
				outfile<<formula("=> " + toBoolVar(s) + " " + adjEdges);
				outfile<<formula("=> " + adjEdges + " " + toBoolVar(s));
			}

			outfile.close();

			time_t t2 = clock();
			if(verbose)
				cout<<"Degree condition assumed"<<endl<<"Time: " + toString(t2 - t1)<<endl;
}

//void Prismatoid::assert_degree()
//{
//	time_t t1 = clock();
//
//	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
//	outfile<<";\tDEGREE CONDITION\n";
//
//	vector<Set> s(nGeo);
//	for(int k = 1; k < dim + 2; k++)
//		for(s[0] = Set::first_n(k); !s[0].contains(nFacet[0]); s[0] = s[0].colex_next())
//			for(s[1] = Set::first_n(dim + 2 - k); !s[1].contains(nFacet[1]); s[1] = s[1].colex_next())
//			{
//				string adjEdges = "(= (+ ";
//				for(int j = 0; j < 2; j++)
//				{
//					if(s[j].cardinality() == 1) continue;
//					for(int i = 0; i < nFacet[j]; i++)
//					{
//						if(!s[j].contains(i)) continue;
//						vector<Set> tmp = s;
//						tmp[j] = s[j] - i;
//						adjEdges += (" " + toVariable(tmp));
//					}
//				}
//				adjEdges += ") ";
//				adjEdges += toString(dim + 1);
//				adjEdges += ")";
//				outfile<<formula("=> " + toBoolVar(s) + " " + adjEdges);
//				//outfile<<formula("=> " + adjEdges + " " + toBoolVar(s));
//			}
//
//			outfile.close();
//
//			time_t t2 = clock();
//			if(verbose)
//				cout<<"Degree condition assumed"<<endl<<"Time: " + toString(t2 - t1)<<endl;
//}
//
void Prismatoid::assert_vertex_on_facet(int geo)
{
	time_t t1 = clock();
	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	outfile<<";\tFACETS HAVE AT LEAST d+1 VERTICES\n";

	for(int f = 0; f < nFacet[geo]; f++)
	{
		int max_sz_bit = 0;
		for(Set s = Set::first_n(dim); !s.contains(nFacet[geo]); s = s.colex_next())
		{
			if(s.contains(f)) continue;

			for(int f2 = 0; f2 < nFacet[1 - geo]; f2++)
				max_sz_bit++;
		}
		max_sz_bit = ceil(log((double)max_sz_bit)/log(2.)); // number of variables in sum

		string extend = "zero_extend[" + toString(max_sz_bit - 1) + "]";
		string bv_ext = "bv" + toString(dim + 1) + "[" + toString(max_sz_bit) + "]";
		string cond = (bit_vector ? "bvuge" : ">=") + string(" (");
		cond += (bit_vector ? "bvadd" : "+");
		for(Set s = Set::first_n(dim); !s.contains(nFacet[geo]); s = s.colex_next())
		{
			if(s.contains(f)) continue;

			for(int f2 = 0; f2 < nFacet[1 - geo]; f2++)
			{
				cond += " ";
				vector<Set> tmp(2);
				tmp[geo] = s + f;
				tmp[1-geo] = Set::singleton(f2);
				cond += "(";
				if(bit_vector)
					cond += extend;
				cond += " ";
				cond += toVariable(tmp) + ")";
			}
		}
		cond += ") ";
		cond += bv_ext;
		outfile<<formula(cond);
	}

	time_t t2 = clock();
	if(verbose)
		cout<<"Vertex on facet assumed"<<endl<<"Time: " + toString(t2 - t1)<<endl;
}

// Works only with d=2 for now
void Prismatoid::assert_vertex_on_face(int d)
{
	time_t t1 = clock();

	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	outfile<<";\tVERTEX ON "<<d<<"-FACE CONDITION\n";

	int max_sz_bit = ceil(log((double)(nFacet[0] + nFacet[1] - dim + 1))/log(2.)); // number of variables in sum
	string extend = "zero_extend[" + toString(max_sz_bit - 1) + "]";
	string bv0_ext = bit_vector ? "bv0[" + toString(max_sz_bit) + "]" : "0";
	string bv2_ext = bit_vector ? "bv2[" + toString(max_sz_bit) + "]" : "2";

	vector<Set> s(nGeo);
	for(int k = 1; k < dim + 1; k++)
		for(s[0] = Set::first_n(k); !s[0].contains(nFacet[0]); s[0] = s[0].colex_next())
			for(s[1] = Set::first_n(dim + 1 - k); !s[1].contains(nFacet[1]); s[1] = s[1].colex_next())
			{
				string sum = bit_vector ? "(bvadd" : "(+";
				for(int i = 0; i < nGeo; i++)
				{
					for(int addToS = 0; addToS < nFacet[i]; addToS++)
					{
						if(s[i].contains(addToS)) continue;
						vector<Set> tmp = s;
						tmp[i] = tmp[i] + addToS;
						sum += " (";
						if(bit_vector)
							sum += extend;
						sum += " " + toVariable(tmp) + ")";
					}
				}
				sum += ")";

				outfile<<formula("=> (= " + sum + " " + bv0_ext + ") " + toBoolVar(s, false));
				outfile<<formula("=> (= " + sum + " " + bv2_ext + ") " + toBoolVar(s));
				outfile<<formula("=> " +  toBoolVar(s, false) + " (= " + sum + " " + bv0_ext + ")");
				outfile<<formula("=> " +  toBoolVar(s) + " (= " + sum + " " + bv2_ext + ")");
			}

			outfile.close();

			time_t t2 = clock();
			if(verbose)
				cout<<"Edge condition assumed"<<endl<<"Time: " + toString(t2 - t1)<<endl;
}

void Prismatoid::get_output()
{
	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	outfile<<":formula (true)\n";
	outfile<<")\n";
	outfile.close();
}

void Prismatoid::assert_number_vertices(int geo, int d, int n)
{
	string str = "= (+ ";
	vector<Set> s(2);
	for(s[geo] = Set::first_n(dim + 1 - d); !s[geo].contains(nFacet[geo]); s[geo] = s[geo].colex_next())
		for(s[1-geo] = Set::first_n(1); !s[1-geo].contains(nFacet[1-geo]); s[1-geo] = s[1-geo].colex_next())
			str += toVariable(s) + " ";

	str += ") " + toString(n);
	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 
	outfile<<formula(str);
}

void Prismatoid::assert_cube(int geo, bool only_bases)
{
	ofstream outfile(file_name.c_str(), std::ios::out | std::ios::binary | std::ios::app); 

	vector<string> v;
	v.push_back("F012x0");
	v.push_back("F023x0");
	v.push_back("F034x1");
	v.push_back("F014x1");

	v.push_back("F125x2");
	v.push_back("F235x2");
	v.push_back("F345x5");
	v.push_back("F145x5");

	v.push_back("F1x024");
	v.push_back("F1x014");
	v.push_back("F1x145");
	v.push_back("F1x245");

	v.push_back("F3x023");
	v.push_back("F3x013");
	v.push_back("F3x135");
	v.push_back("F3x235");

	if(!only_bases)
	{
		v.push_back("F34x15");
		v.push_back("F14x15");
		v.push_back("F01x01");
		v.push_back("F03x01");

		v.push_back("F12x02");
		v.push_back("F23x02");
		v.push_back("F35x25");
		v.push_back("F15x25");
	}
	string str;
	for(int i = 0; i < v.size(); i++)
		str += v[i] + "|";
	//outfile<<str<<"\n";

	if(!only_bases)
	{
		vector<Set> s(2);
		for(int k = 1; k < dim + 2; k++)
		{
			for(s[0] = Set::first_n(k); !s[0].contains(nFacet[0]); s[0] = s[0].colex_next())
				for(s[1] = Set::first_n(dim + 2 - k); !s[1].contains(nFacet[1]); s[1] = s[1].colex_next())
				{
					if(find(v.begin(), v.end(), toVariable(s)) != v.end())
						outfile<<formula(toBoolVar(s, true));
					else
						outfile<<formula(toBoolVar(s, false));
				}	
		}
	}
	else
	{
		for(int i = 0; i < 2; i++)
		{
			vector<Set> s(2);
			for(s[i] = Set::first_n(1); !s[i].contains(nFacet[i]); s[i] = s[i].colex_next())
				for(s[1-i] = Set::first_n(dim + 2 - 1); !s[1-i].contains(nFacet[1-i]); s[1-i] = s[1-i].colex_next())
				{
					if(find(v.begin(), v.end(), toVariable(s)) != v.end())
						outfile<<formula(toBoolVar(s, true));
					else
						outfile<<formula(toBoolVar(s, false));
				}	
		}
	}
}
