#include "Utils/Types.h"
#include "Boundary/BoundaryCondition.h"


void parse(std::istream& stream,
               List<Vector>& vertices,
               List<Array<Index,2>>& faces,
               List<List<Index>>& cells,
               HashMap<Index, Boundaries>& boundaries);


List<std::string> split(std::string const& s, char sep = ' ');