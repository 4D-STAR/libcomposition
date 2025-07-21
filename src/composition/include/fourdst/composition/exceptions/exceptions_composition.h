#pragma once

#include <exception>
#include <string>
#include <iostream>

namespace fourdst::composition::exceptions {
    class CompositionError : public std::exception {
    protected:
        std::string m_message;
    public:
        explicit CompositionError(const std::string& message)
            : m_message(std::move(message)) {}

        const char* what() const noexcept override{
            return m_message.c_str();
        }
    };

    class CompositionEntryError : public std::exception {
    protected:
        std::string m_message;
    public:
        explicit CompositionEntryError(const std::string& message)
            : m_message(std::move(message)) {}

        const char* what() const noexcept override {
            return m_message.c_str();
        }
    };

    class CompositionNotFinalizedError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    class InvalidCompositionError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    class InvalidMixingMode final : public CompositionError {
        using CompositionError::CompositionError;
    };

    class InvalidSymbolError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    class UnregisteredSymbolError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    class FailedToFinalizeCompositionError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    class InvalidSpeciesSymbolError final : public CompositionEntryError {
        using CompositionEntryError::CompositionEntryError;
    };

    class EntryAlreadyInitializedError final : public CompositionEntryError {
        using CompositionEntryError::CompositionEntryError;
    };

    class CompositionModeError final : public CompositionEntryError {
        using CompositionEntryError::CompositionEntryError;
    };

}