#include "Parse.h"


List<std::string> split(std::string const& s, char sep)
{
    List<std::string> res;

    size_t i = 0;
    while (i < s.size())
    {
        while (i < s.size() && s[i] == sep)
            i++;
        size_t j = i;
        while (j < s.size() && s[j] != sep)
            j++;
        if (j > i)
            res.push_back(s.substr(i, j-i));
        i = j;
    }
    
    return res;
}


static void throwParseError(Index lineIdx);
static void parseVertex(List<std::string> const& tokens, List<Vector>& vertices);
static void parseFace(List<std::string> const& tokens, List<Array<Index,2>>& faces);
static void parseCell(List<std::string> const& tokens, List<List<Index>>& cells);
static bool parseBoundaries(List<std::string> const& tokens, HashMap<Index, Boundaries>& boundariesMap);


void parse(std::istream& stream,
               List<Vector>& vertices,
               List<Array<Index,2>>& faces,
               List<List<Index>>& cells,
               HashMap<Index, Boundaries>& boundariesMap)
{
    Index lineIdx = 0;

    for (; !stream.eof(); lineIdx++)
    {
        std::string line;
        std::getline(stream, line);

        List<std::string> tokens = split(line, ' ');

        if (tokens.size() == 0 || std::string("#") == tokens[0])
            continue;
        
        if (tokens[0].size() > 1)
            throwParseError(lineIdx);

        char entity = tokens[0].at(0);
        switch (entity)
        {
        case 'v':
        {
            parseVertex(tokens, vertices);
            break;
        }
        case 'f':
        {
            parseFace(tokens, faces);
            break;
        }
        case 'c':
        {
            parseCell(tokens, cells);
            break;
        }
        case 'b':
        {
            if (!parseBoundaries(tokens, boundariesMap))
                throwParseError(lineIdx);
            break;
        }
        default:
            throwParseError(lineIdx);
        }       
    }
}


static void parseVertex(List<std::string> const& tokens, List<Vector>& vertices)
{
    Scalar x = std::stod(tokens.at(1));
    Scalar y = std::stod(tokens.at(2));
    vertices.push_back({x,y,0});
}


static void parseFace(List<std::string> const& tokens, List<Array<Index,2>>& faces)
{
    Index v1 = std::stoi(tokens.at(1));
    Index v2 = std::stoi(tokens.at(2));
    faces.push_back({v1, v2});
}


static void parseCell(List<std::string> const& tokens, List<List<Index>>& cells)
{
    cells.push_back({});
    for (Index idx = 1; idx < (Index)tokens.size(); idx++)
        cells.back().push_back(std::stoi(tokens[idx]));
}


static bool parseBoundaries(List<std::string> const& tokens, HashMap<Index, Boundaries>& boundariesMap)
{
    Index faceIdx = std::stoi(tokens.at(1));
    Boundaries& boundaries = boundariesMap[faceIdx];

    using enum BoundaryConditionType;
    BoundaryConditionType boundariesType;

    if (tokens.at(3) == std::string("fixedValue"))
        boundariesType = FIXED_VALUE;
    else if (tokens.at(3) == std::string("fixedGradient"))
        boundariesType = FIXED_GRADIENT;
    else
        return false;

    if (tokens.at(2) == std::string("U"))
    {
        boundaries.uBoundary.type = boundariesType;
        boundaries.uBoundary.value = {std::stod(tokens.at(4)), std::stod(tokens.at(5)), 0};
    }
    else if (tokens.at(2) == std::string("p"))
    {
        boundaries.pBoundary.type = boundariesType;
        boundaries.pBoundary.value = std::stod(tokens.at(4));
    }
    else
        return false;
    
    return true;
}


static void throwParseError(Index lineIdx)
{
    throw std::runtime_error("Unrecognized symbol while parsing PolyMesh2D, line " + std::to_string(lineIdx+1));
}