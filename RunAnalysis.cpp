#include "LLGAnalysis.h"
#include <vector>
#include <string>

int main( int argc, char **argv ) {
    LLGAnalysis *analysis = LLGAnalysis::GetInstance( argv[1] );
    cout << "now initing" << endl;
    analysis->Init();
    cout << "starting evt loop" << endl;
    analysis->RunEventLoop( (argc > 2 ) ? atoi(argv[2]) : -1 );
    cout << "\nfinishing " << endl;
    analysis->FinishRun();
    return 0;
}
