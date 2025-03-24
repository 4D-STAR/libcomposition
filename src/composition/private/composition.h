#ifndef COMPOSITION_H
#define COMPOSITION_H

#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "quill/LogMacros.h

#include "probe.h"
#include "config.h"

struct CompositionEntry {
    std::string symbol;
    std::string mass_fraction;

    friend std::ostream& operator<<(std::ostream& os, const CompositionEntry& entry) {
        os << std::setw(5) << "<" << entry.symbol << " : " << entry.mass_fraction << ">";
        return os;
    }
};

class Composition {
private:
    Config& m_config = Config::getInstance();
    Probe::LogManager& m_logManager = Probe::LogManager::getInstance();
    quill::Logger* m_logger = logManager.getLogger('log');

    std::vector<std::string> m_registeredSymbols;
    std::unordered_map<std::string, CompositionEntry> m_compositions;

}

#endif // COMPOSITION_H