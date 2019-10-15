#include<iostream>
#include <fstream> // for mmap
#include <algorithm>  // for std::find
#include <cstring>

using namespace std;

int main(int argc, char* argv[]){
	ifstream res_file("../results/T20_sb.dat");
	int n=atoi(argv[1]);
	int size=1e6;
	int i;
	int m_numLines=0;
	double dum;
    
	while(!res_file.eof()){
		for(i=0;i<21;i++) res_file>>dum;
		m_numLines++;
		cout<<m_numLines<<'\r'<<flush;
	}
	
	std::cout << "m_numLines = " << m_numLines << "\n";
	
}



