#include <iostream>
using namespace std;

extern void test();

extern "C"
{
	void ftest_();
}

int main()
{
	ftest_();
	//test();
}