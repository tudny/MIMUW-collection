#include "indexer/indexer.h"
#include "query/query.h"
#include "logger.h"
#include "params.h"
#include "error.h"

int main(int argc, char *argv[]) {

    LOGDEBUG("Running the main function");

    InputParams params = parse_args(argc, argv);

    std::visit(overloaded{
            [](const IndexingParams &params) { index(params); },
            [](const QueryParams &params) { query(params); },
            [](std::monostate) { print_usage(); }
    }, params);

    return 0;
}
