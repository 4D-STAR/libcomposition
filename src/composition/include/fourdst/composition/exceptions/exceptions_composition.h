#pragma once

#include <exception>
#include <string>

namespace fourdst::composition::exceptions {
    /**
     * @class CompositionError
     * @brief Base class for exceptions related to composition objects.
     *
     * This exception is thrown when an error occurs at the composition level,
     * such as invalid configuration or state.
     */
    class CompositionError : public std::exception {
    protected:
        /**
         * @brief The error message.
         */
        std::string m_message;
    public:
        /**
         * @brief Constructs a CompositionError with an error message.
         * @param message The error message.
         */
        explicit CompositionError(const std::string& message)
            : m_message(std::move(message)) {}

        /**
         * @brief Returns the error message.
         * @return A C-style string containing the error message.
         */
        const char* what() const noexcept override{
            return m_message.c_str();
        }
    };

    /**
     * @class CompositionEntryError
     * @brief Base class for exceptions related to individual entries within a composition.
     *
     * This exception is thrown for errors specific to a single component or entry
     * in a composition, such as an invalid species symbol or duplicate initialization.
     */
    class CompositionEntryError : public std::exception {
    protected:
        /**
         * @brief The error message.
         */
        std::string m_message;
    public:
        /**
         * @brief Constructs a CompositionEntryError with an error message.
         * @param message The error message.
         */
        explicit CompositionEntryError(const std::string& message)
            : m_message(std::move(message)) {}

        /**
         * @brief Returns the error message.
         * @return A C-style string containing the error message.
         */
        const char* what() const noexcept override {
            return m_message.c_str();
        }
    };

    /**
     * @class CompositionNotFinalizedError
     * @brief Exception thrown when an operation is attempted on a composition that has not been finalized.
     *
     * Certain operations require the composition to be in a "finalized" state.
     * This error indicates that such an operation was called prematurely.
     */
    class CompositionNotFinalizedError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    /**
     * @class InvalidCompositionError
     * @brief Exception thrown when a composition is in an invalid or inconsistent state.
     */
    class InvalidCompositionError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    /**
     * @class InvalidMixingMode
     * @brief Exception thrown for an invalid or unsupported mixing mode.
     *
     * Compositions can be defined with different mixing modes (e.g., by mass, by mole).
     * This error is thrown if an invalid mode is specified.
     */
    class InvalidMixingMode final : public CompositionError {
        using CompositionError::CompositionError;
    };

    /**
     * @class InvalidSymbolError
     * @brief Exception thrown when a symbol used in a composition is invalid.
     */
    class InvalidSymbolError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    /**
     * @class UnregisteredSymbolError
     * @brief Exception thrown when a symbol is used that has not been registered.
     *
     * This typically occurs when a chemical species is used that is not known to the system.
     */
    class UnregisteredSymbolError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    /**
     * @class FailedToFinalizeCompositionError
     * @brief Exception thrown when the finalization process of a composition fails.
     */
    class FailedToFinalizeCompositionError final : public CompositionError {
        using CompositionError::CompositionError;
    };

    /**
     * @class InvalidSpeciesSymbolError
     * @brief Exception thrown for an invalid chemical species symbol in a composition entry.
     */
    class InvalidSpeciesSymbolError final : public CompositionEntryError {
        using CompositionEntryError::CompositionEntryError;
    };

    /**
     * @class EntryAlreadyInitializedError
     * @brief Exception thrown when attempting to initialize a composition entry that has already been initialized.
     */
    class EntryAlreadyInitializedError final : public CompositionEntryError {
        using CompositionEntryError::CompositionEntryError;
    };

    /**
     * @class CompositionModeError
     * @brief Exception thrown due to a conflict in composition modes at the entry level.
     *
     * This may occur if an entry's configuration is incompatible with the overall composition's mode.
     */
    class CompositionModeError final : public CompositionEntryError {
        using CompositionEntryError::CompositionEntryError;
    };

}