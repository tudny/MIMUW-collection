#include <iostream>
#include <regex>
#include <set>
#include <queue>
#include <unordered_map>
#include <unordered_set>

using namespace std;

namespace {
    namespace gates {
        enum gate_t {
            NOT,
            XOR,
            AND,
            NAND,
            OR,
            NOR,
            NOT_GATE_ERROR
        };

        static const array<gate_t, 6> gateTypes {NOT, XOR, AND, NAND, OR, NOR};

        bool andFunction(bool p, bool q) {
            return p && q;
        }

        bool orFunction(bool p, bool q) {
            return p || q;
        }

        bool xorFunction(bool p, bool q) {
            return p ^ q;
        }

        bool notFunction([[maybe_unused]] bool p, bool q) {
            return !q;
        }

        using gate_function_t = bool (*)(bool, bool);
        // element neutralny, negacja wyniku funkcji, przerywalna
        using gate_group_t = tuple<gate_function_t, bool, bool, bool>;
        using gate_functions_map_t = unordered_map<gate_t, gate_group_t>;

        gate_functions_map_t initializeGateFunctionsMap() {
            gate_functions_map_t gateFunctions;

            gateFunctions[AND] = make_tuple(andFunction, true, false, true);
            gateFunctions[OR] = make_tuple(orFunction, false, false, true);
            gateFunctions[NAND] = make_tuple(andFunction, true, true, true);
            gateFunctions[NOR] = make_tuple(orFunction, false, true, true);
            gateFunctions[XOR] = make_tuple(xorFunction, false, false, false);
            gateFunctions[NOT] = make_tuple(notFunction, false, false, false);

            return gateFunctions;
        }

        gate_group_t &getGroupOfGate(const gate_t &gate) {
            static gate_functions_map_t gateFunctions = initializeGateFunctionsMap();
            return gateFunctions[gate];
        }

        using gate_to_regex_string_t = unordered_map<gate_t, string>;
        using gate_to_regex_t = unordered_map<gate_t, regex>;

        gate_to_regex_string_t initializeGateToRegexMap() {
            gate_to_regex_string_t gateRegexStringMap;

            const static string numberRegex = R"(([1-9]\d{0,8}))";

            unordered_map<gate_t, pair<string, string>> regexBuilderBlocks;
            regexBuilderBlocks[NOT] = {"NOT", "2"};
            regexBuilderBlocks[XOR] = {"XOR", "3"};
            regexBuilderBlocks[AND] = {"AND", "3,"};
            regexBuilderBlocks[NAND] = {"NAND", "3,"};
            regexBuilderBlocks[OR] = {"OR", "3,"};
            regexBuilderBlocks[NOR] = {"NOR", "3,"};

            for (const auto &[gate, data] : regexBuilderBlocks) {
                string regexStr = "^\\s*" + data.first + "(\\s+" + numberRegex + "){" + data.second + "}\\s*$";
                gateRegexStringMap[gate] = regexStr;
            }

            return gateRegexStringMap;
        }

        regex createRegexOfGate(const gate_t &gate) {
            static unordered_map<gate_t, string> gateRegexStringMap = initializeGateToRegexMap();

            regex gateRegex(gateRegexStringMap[gate]);
            return gateRegex;
        }

        gate_to_regex_t initializeRegexes() {
            gate_to_regex_t gateRegexMap;

            for (const gate_t &gate : gateTypes)
                gateRegexMap[gate] = createRegexOfGate(gate);

            return gateRegexMap;
        }

        gate_t getGateTypeOfString(const string &str) {
            static gate_to_regex_t gateRegexMap = initializeRegexes();

            auto it = find_if(gateTypes.begin(), gateTypes.end(), [&str](const gate_t &gateType) {
                const regex gateRegex = gateRegexMap.at(gateType);
                return regex_match(str, gateRegex);
            });

            return it == gateTypes.end() ? NOT_GATE_ERROR : *it;
        }
    }

    namespace validation {
        void removeDuplicateSpaces(string &str) {
            auto strEnd = unique(str.begin(), str.end(), [](char prev, char next) {
                return isspace(prev) && isspace(next);
            });

            str.erase(strEnd, str.end());
        }

        void whitespacesToSpaces(string &str) {
            replace_if(str.begin(), str.end(), [](char c) {
                return isspace(c);
            }, ' ');
        }

        void trimLeft(string &str) {
            auto new_begin = find_if(str.begin(), str.end(), [](char c) {
                return !isspace(c);
            });

            str.erase(str.begin(), new_begin);
        }

        void trimRight(string &str) {
            auto new_end = find_if(str.rbegin(), str.rend(), [](char c) {
                return !isspace(c);
            });

            str.erase(new_end.base(), str.end());
        }

        void trimSpaces(string &str) {
            trimRight(str);
            trimLeft(str);
        }

        void cleanString(string &str) {
            removeDuplicateSpaces(str);
            whitespacesToSpaces(str);
            trimSpaces(str);
        }

        bool isProperLine(const string &cleanStr, gates::gate_t &gate) {
            return (gate = gates::getGateTypeOfString(cleanStr)) != gates::NOT_GATE_ERROR;
        }
    }

    namespace graph {
        using node_t = uint32_t;

        struct node_hasher_t {
            size_t operator()(const node_t &node) const {
                static node_t magic = 20211014;
                return node ^ magic;
            }
        };

        using node_list_t = set<node_t>;
        using node_inputs_t = unordered_multiset<node_t, node_hasher_t>;
        using node_list_desc_t = set<node_t, greater<>>;
        using node_value_map_t = map<node_t, bool>;
        using node_meta_t = tuple<node_inputs_t, gates::gate_t>;
        using gate_nodes_t = unordered_map<node_t, node_meta_t, node_hasher_t>;
        using graph_t = pair<node_list_t, gate_nodes_t>;

        node_list_t &getAllNodes(graph_t &graph) {
            return get<0>(graph);
        }

        const node_list_t &getAllNodes(const graph_t &graph) {
            return get<0>(graph);
        }

        gate_nodes_t &getGateNodes(graph_t &graph) {
            return get<1>(graph);
        }

        const gate_nodes_t &getGateNodes(const graph_t &graph) {
            return get<1>(graph);
        }

        node_inputs_t &getNeighbours(node_meta_t &nodeData) {
            return get<0>(nodeData);
        }

        const node_inputs_t &getNeighbours(const node_meta_t &nodeData) {
            return get<0>(nodeData);
        }

        const gates::gate_t &getGate(const node_meta_t &nodeData) {
            return get<1>(nodeData);
        }

        void addNode(graph_t &graph, const node_t &node) {
            getAllNodes(graph).insert(node);
        }

        bool addGateNode(graph_t &graph, const node_t node, const gates::gate_t gate) {
            if (graph::getGateNodes(graph).contains(node))
                return false;

            node_inputs_t emptyList;
            graph::getGateNodes(graph).insert({node, make_tuple(emptyList, gate)});

            addNode(graph, node);

            return true;
        }

        void addEdge(graph_t &graph, const node_t begin, const node_t end) {
            addNode(graph, begin);
            getNeighbours(graph::getGateNodes(graph).at(end)).insert(begin);
        }

        namespace algorithm {

            using node_counter_t = unordered_map<node_t, size_t>;
            using node_order_t = vector<node_t>;
            using node_queue_t = queue<node_t>;

            void topoSortPreprocessing(const graph_t &graph, node_counter_t &howManyInputs,
                                       node_queue_t &availableNodes) {

                auto &allNodes = graph::getAllNodes(graph);

                for (const auto &node : allNodes) {
                    howManyInputs.insert({node, 0});
                }

                for (auto &[_, data] : graph::getGateNodes(graph)) {
                    for (const auto &neighbor : getNeighbours(data)) {
                        howManyInputs[neighbor]++;
                    }
                }

                for (const auto &node: allNodes) {
                    if (!howManyInputs[node]) {
                        availableNodes.push(node);
                    }
                }
            }

            void processAvailableNodes(const graph_t &graph, node_queue_t &availableNodes,
                                       node_counter_t &howManyInputs,
                                       node_order_t &toposortedNodes) {

                while (!availableNodes.empty()) {
                    const node_t node = availableNodes.front();
                    availableNodes.pop();

                    toposortedNodes.push_back(node);

                    auto &gateNodes = graph::getGateNodes(graph);
                    if (gateNodes.contains(node)) {
                        for (const auto &nextNode : getNeighbours(gateNodes.at(node))) {
                            if (--howManyInputs[nextNode] == 0) {
                                availableNodes.push(nextNode);
                            }
                        }
                    }
                }
            }

            node_order_t topoSort(graph_t &graph) {
                node_counter_t howManyInputs;
                node_queue_t availableNodes;
                node_order_t toposortedNodes;

                topoSortPreprocessing(graph, howManyInputs, availableNodes);
                processAvailableNodes(graph, availableNodes, howManyInputs, toposortedNodes);

                reverse(toposortedNodes.begin(), toposortedNodes.end());

                return toposortedNodes;
            }

            bool hasCycle(const graph_t &graph, const node_order_t &topologicalOrder) {
                return graph::getAllNodes(graph).size() != topologicalOrder.size();
            }
        }
    }

    namespace simulation {
        using bit_mask_t = u_int64_t;

        void nextBitMask(bit_mask_t &bitMask) {
            bitMask++;
        }

        bool isBitOn(const bit_mask_t bitMask, const size_t n) {
            return ((1 << n) & bitMask) != 0;
        }

        void printResult(const graph::graph_t &graph, const graph::node_value_map_t &allNodesValues) {
            for (const auto &node : graph::getAllNodes(graph)) {
                cout << allNodesValues.at(node);
            }
            cout << "\n";
        }

        graph::node_list_desc_t getInputSignals(const graph::graph_t &graph) {
            graph::node_list_desc_t inputSignals;

            for (const auto &node : graph::getAllNodes(graph)) {
                if (!graph::getGateNodes(graph).contains(node)) {
                    inputSignals.insert(node);
                }
            }

            return inputSignals;
        }

        void setInputSignalsValues(const graph::node_list_desc_t &inputSignals, const bit_mask_t &inputValues,
                                   graph::node_value_map_t &allNodesValues) {
            size_t index = 0;
            for (const auto &node : inputSignals) {
                allNodesValues[node] = isBitOn(inputValues, index++);
            }
        }

        bool calculateNodeValue(const graph::graph_t &graph, const graph::node_t &node,
                                const graph::node_value_map_t &allNodesValues) {

            auto &metaData = graph::getGateNodes(graph).at(node);
            auto &neighbours = graph::getNeighbours(metaData);
            auto gateType = graph::getGate(metaData);
            auto &[function, neutralElement, shouldNegate, breakable] = gates::getGroupOfGate(gateType);

            bool nodeValue = neutralElement;
            for (const auto &neighbourNode : neighbours) {
                nodeValue = function(nodeValue, allNodesValues.at(neighbourNode));
                if (breakable && nodeValue != neutralElement) break;
            }

            return shouldNegate == !nodeValue;
        }

        void runSingleLayout(const graph::graph_t &graph, const graph::node_list_desc_t &inputSignals,
                             const bit_mask_t &inputValues, const graph::algorithm::node_order_t &toposorted) {

            graph::node_value_map_t allNodesValues;

            setInputSignalsValues(inputSignals, inputValues, allNodesValues);

            for (const auto &node : toposorted) {
                if (!inputSignals.contains(node)) {
                    allNodesValues[node] = calculateNodeValue(graph, node, allNodesValues);
                }
            }

            printResult(graph, allNodesValues);
        }

        void run(graph::graph_t &graph) {
            auto toposorted = graph::algorithm::topoSort(graph);

            if (graph::algorithm::hasCycle(graph, toposorted)) {
                cerr << "Error: sequential logic analysis has not yet been implemented.\n";
                return;
            }

            auto inputSignals = getInputSignals(graph);
            gates::initializeGateFunctionsMap();

            bit_mask_t inputValues = 0;
            const bit_mask_t allInputsOn = (1 << inputSignals.size()) - 1;

            while (inputValues <= allInputsOn) {
                runSingleLayout(graph, inputSignals, inputValues, toposorted);
                nextBitMask(inputValues);
            }
        }
    }

    namespace input {
        void printError(size_t lineNumber, string &message) {
            cerr << "Error in line " << lineNumber << ": " << message << "\n";
        }

        /**
         * @return jeśli jest error, to numer node'a, który powoduje zwracie. wpp 0.
         */
        graph::node_t analizeInput(const string &input, gates::gate_t &gate, graph::graph_t &graph) {
            stringstream stream(input);
            string str;

            stream.ignore(numeric_limits<streamsize>::max(), ' ');
            stream >> str;
            graph::node_t outputId = stoi(str);

            if (graph::addGateNode(graph, outputId, gate)) {
                while (stream) {
                    stream >> str;
                    if (!str.empty()) {
                        graph::node_t inputId = stoi(str);
                        graph::addEdge(graph, inputId, outputId);
                        str.clear();
                    }
                }
                return 0;
            } else
                return outputId;
        }

        bool readInput(graph::graph_t &graph) {
            string input;
            size_t lineNumber = 1;
            bool errorDetected = false;

            while (getline(cin, input)) {
                string originalInput = input;
                validation::cleanString(input);

                gates::gate_t gate;

                if (validation::isProperLine(input, gate)) {
                    graph::node_t outputId = analizeInput(input, gate, graph);
                    if (outputId != 0) {
                        string message = "signal " + to_string(outputId) +
                                         " is assigned to multiple outputs.";
                        errorDetected = true;
                        printError(lineNumber, message);
                    }
                } else {
                    errorDetected = true;
                    printError(lineNumber, originalInput);
                }

                lineNumber++;
            }

            return !errorDetected;
        }
    }
}

int main() {
    graph::graph_t graph;
    if (input::readInput(graph)) {
        simulation::run(graph);
    }

    return 0;
}
