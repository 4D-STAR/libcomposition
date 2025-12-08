#pragma once

#include <exception>
#include <string>
#include <utility>

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
        explicit CompositionError(std::string  message)
            : m_message(std::move(message)) {}

        /**
         * @brief Returns the error message.
         * @return A C-style string containing the error message.
         */
        [[nodiscard]] const char* what() const noexcept override{
            return m_message.c_str();
        }
    };

    /**
     * @class InvalidCompositionError
     * @brief Exception thrown when a composition is in an invalid or inconsistent state.
     */
    class InvalidCompositionError final : public CompositionError {
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
     * @class SpeciesError
     * @brief Base class for exceptions related to atomic species.
     */
    class SpeciesError : public std::exception {
    protected:
        std::string m_message;
    public:
        explicit SpeciesError(std::string  message)
            : m_message(std::move(message)) {}

        [[nodiscard]] const char* what() const noexcept override {
            return m_message.c_str();
        }
    };

    /**
     * @class UnknownSymbolError
     * @brief Exception thrown when an unknown symbol is encountered.
     *
     * This typically occurs when a symbol does not correspond to any known atomic species.
     */
    class UnknownSymbolError final : public SpeciesError {
        using SpeciesError::SpeciesError;
    };

}