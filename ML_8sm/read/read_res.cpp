#include<iostream>
#include <boost/iostreams/device/mapped_file.hpp> // for mmap
#include <algorithm>  // for std::find
#include <cstring>

using namespace std;

int main(int argc, char* argv[]){
	//ifstream res_file("../results/T20_sb.dat");
	int n=atoi(argv[1]);
	int size=1e6;
	int i;

    boost::iostreams::mapped_file mmap("../results/T20_sb.dat", boost::iostreams::mapped_file::readonly);
	auto f = mmap.const_data();
	auto l = f + mmap.size();
	
	uintmax_t m_numLines = 0;
	while (f && f!=l)
		if ((f = static_cast<const char*>(memchr(f, '\n', l-f))))
			m_numLines++, f++;
			
		std::cout << "m_numLines = " << m_numLines << "\n";
	
}



