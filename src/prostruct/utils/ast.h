/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_AST_H
#define PROSTRUCT_AST_H

#define FMT_STRING_ALIAS 1

#include <armadillo>

#include <fmt/format.h>

namespace prostruct
{
	template <typename T>
	class StructBase;

	template <typename T>
	class Atom;

	namespace parser
	{
		namespace detail
		{
			/**
			 * The exception raised by any class involved
			 * in the parsing and interpreting of
			 * Prostruct's DSL
			 */
			class ParserException : public std::exception
			{
			public:
				/**
				 * The ParserException constructor that can
				 * be passed a messaged and an arbitrary number
				 * of parameters.
				 *
				 * @tparam Args the argument pack types
				 * @param message the main message that is passed to
				 * fmt::format
				 * @param args the additional values passed to the
				 * main message
				 */
				template <typename... Args>
				ParserException(const std::string& message, Args... args)
				{
					m_message = fmt::format(message, args...);
				}
				/**
				 * Returns the message as a char array
				 *
				 * @return the exception message
				 */
				const char* what() const noexcept final
				{
					return format(fmt("{}{}"), m_error_msg_start, m_message).c_str();
				}

				/** The start of the error message */
				static constexpr std::string_view m_error_msg_start = "Parser error: ";

			private:
				/** The exception message */
				std::string m_message;
			};

			/**
			 * The possible token types. Note that this enum
			 * is also used by the lexer to merger ALPHA and
			 * NUMERIC to form ATOM_NAME, ATOM_NUMBER and so
			 * on
			 */
			enum class TOKEN_TYPE
			{
				NONE,
				NUMERIC,
				ALPHA,
				RPAREN,
				LPAREN,
				AND,
				OR,
				RANGE,
				ATOM_NAME,
				ATOM_NUMBER,
				RESIDUE_NAME,
				RESIDUE_NUMBER,
				CHAIN_NAME,
				EOF_
			};

			/**
			 * A helper function to make TOKEN_TYPE human readable
			 * @param token_type the token type enum
			 * @return the token type string representation
			 */
			inline std::string get_string_from_enum(TOKEN_TYPE token_type) noexcept
			{
				switch (token_type)
				{
				case TOKEN_TYPE::NONE:
				{
					return "NONE";
				}
				case TOKEN_TYPE::NUMERIC:
				{
					return "NUMERIC";
				}
				case TOKEN_TYPE::ALPHA:
				{
					return "ALPHA";
				}
				case TOKEN_TYPE::RPAREN:
				{
					return "RIGHT_PARENTHESIS";
				}
				case TOKEN_TYPE::LPAREN:
				{
					return "LEFT_PARENTHESIS";
				}
				case TOKEN_TYPE::ATOM_NAME:
				{
					return "KEYWORD_ATOM_NAME";
				}
				case TOKEN_TYPE::RESIDUE_NAME:
				{
					return "KEYWORD_RESIDUE_NAME";
				}
				case TOKEN_TYPE::RESIDUE_NUMBER:
				{
					return "KEYWORD_RESIDUE_NUMBER";
				}
				case TOKEN_TYPE::CHAIN_NAME:
				{
					return "KEYWORD_CHAIN_NAME";
				}
				case TOKEN_TYPE::AND:
				{
					return "LOGICAL_AND";
				}
				case TOKEN_TYPE::OR:
				{
					return "LOGICAL_OR";
				}
				case TOKEN_TYPE::RANGE:
				{
					return "RANGE";
				}
				case TOKEN_TYPE::EOF_:
				{
					return "EOF";
				}
				default:
				{
					return "N/A";
				}
				}
			}

			/**
			 * Extension to the enum class TOKEN_TYPE ++operator
			 * This is needed to allow switching from NAME -> NUMBER
			 * with ++TOKEN_TYPE
			 */
			inline TOKEN_TYPE& operator++(TOKEN_TYPE& lhs)
			{
				lhs = static_cast<TOKEN_TYPE>(
					static_cast<std::underlying_type_t<TOKEN_TYPE>>(lhs) + 1);
				return lhs;
			}

			/**
			 * Convert a string to lower case.
			 * Straight from https://en.cppreference.com/w/cpp/string/byte/tolower.
			 * @param s
			 * @return
			 */
			inline std::string str_tolower(std::string s)
			{
				std::transform(
					s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
				return s;
			}

			/**
			 * A token maps a value to its type
			 * and position in the expression
			 */
			class Token
			{
			public:
				/**
				 * The default Token constructor used as a placeholder
				 */
				Token()
					: m_type(TOKEN_TYPE::NONE)
					, m_value("")
					, m_position(0)
				{
				}
				/**
				 * The Token constructor that maps a value to
				 * its type and the records its position in the
				 * expression.
				 *
				 * @param name the token type
				 * @param value the value of the token
				 * @param position the position in the expression
				 */
				Token(TOKEN_TYPE name, const std::string& value, size_t position)
					: m_type(name)
					, m_value(value)
					, m_position(position)
				{
				}

				/** Token type getter */
				TOKEN_TYPE get_type() const noexcept { return m_type; }
				/** Token value getter */
				std::string get_value() const noexcept { return m_value; }
				/** Token position getter */
				size_t get_position() const noexcept { return m_position; }
				/** Print the token representation */
				void print() const noexcept { std::cout << get_repr() << "\n"; }
				/** Returns the Token representation as a string */
				std::string get_repr() const noexcept
				{
					return format(fmt("Token(value={}, type={}, pos={})"), m_value,
						get_string_from_enum(m_type), m_position);
				}

			private:
				/** The Token type */
				TOKEN_TYPE m_type;
				/** The Token value */
				std::string m_value;
				/** The Token position in the expression */
				size_t m_position;
			};

			/**
			 * The Node in the AST tree with various
			 * helper functions.
			 */
			class Node
			{
			public:
				/**
				 * Node default constructor with a
				 * default Token
				 */
				Node() {}
				/**
				 * Node constructor for the tree stumps,
				 * i.e. does not have a left or right Node
				 * connection.
				 *
				 * @param token the Token instance of the node
				 */
				Node(Token token)
					: m_token(std::move(token))
				{
				}
				/**
				 * Node constructor with further Node
				 * connections.
				 *
				 * @param left the left Node
				 * @param token the Token instance of the node
				 * @param right the right Node
				 */
				Node(const std::shared_ptr<Node>& left, Token token,
					const std::shared_ptr<Node>& right)
					: m_left(left)
					, m_token(std::move(token))
					, m_right(right)
				{
					if (m_token.get_type() != TOKEN_TYPE::OR
						&& m_token.get_type() != TOKEN_TYPE::AND)
						throw ParserException(
							"Expected op to be either of type {} or {}, but got {}",
							get_string_from_enum(TOKEN_TYPE::OR),
							get_string_from_enum(TOKEN_TYPE::AND),
							get_string_from_enum(m_token.get_type()));
				}

				/** Returns the operation type given by the Token */
				TOKEN_TYPE get_op_type() const noexcept { return m_token.get_type(); }
				/** Returns the string representation of the right Node */
				std::string get_right_repr() const noexcept
				{
					if (m_right)
						return get_right()->get_repr();
					else
						return std::string("N/A");
				}
				/** Returns the string representation of the left Node */
				std::string get_left_repr() const noexcept
				{
					if (m_right)
						return get_left()->get_repr();
					else
						return std::string("N/A");
				}
				/** Returns the string representation of this Node */
				std::string get_repr() const noexcept
				{
					return format(fmt("Node(left={}, right={}, token={})"), get_left_repr(),
						get_right_repr(), m_token.get_repr());
				}
				/** Left Node getter */
				std::shared_ptr<Node> get_left() const noexcept { return m_left; }
				/** Right Node getter */
				std::shared_ptr<Node> get_right() const noexcept { return m_right; }
				/** Token type getter */
				TOKEN_TYPE get_type() const noexcept { return m_token.get_type(); }
				/** Token value getter */
				std::string get_value() const noexcept { return m_token.get_value(); }

			private:
				/** the left Node */
				std::shared_ptr<Node> m_left;
				/** the right Node */
				std::shared_ptr<Node> m_right;
				/** the node Token */
				Token m_token;
			};

			/**
			 * The Lexer parses a string and splits it into
			 * Tokens. The Lexer ignores whitespaces.
			 */
			class Lexer
			{
			public:
				/**
				 * Lexer constructor for a given
				 * expression.
				 *
				 * @param text the text to be tokenised
				 */
				Lexer(const std::string& text)
					: m_text_beginning(text.cbegin())
					, m_text_start(text.cbegin())
					, m_text_iter(text.cbegin())
					, m_text_end(text.cend())
					, internal_token(TOKEN_TYPE::NONE)
				{
				}

				/**
				 * Returns the next token. First call starts
				 * at the beginning of the sentence.
				 * The end of the expression return a EOF Token.
				 *
				 * @return the next Token instance
				 */
				Token get_next_token()
				{
					m_text_start = m_text_iter;

					// skip whitespace
					while (std::isspace(*m_text_iter))
					{
						++m_text_iter;
						m_text_start = m_text_iter;
					}
					// return EOF
					if (m_text_iter == m_text_end)
					{
						return Token(
							TOKEN_TYPE::EOF_, "EOF", std::distance(m_text_beginning, m_text_start));
					}
					// process punctuation
					if (valid_punctuation())
					{
						if (*m_text_iter == ')')
							internal_token = TOKEN_TYPE::RPAREN;
						else if (*m_text_iter == '(')
							internal_token = TOKEN_TYPE::LPAREN;
						else
							throw ParserException("This is a bug. A new punctuation has been added "
												  "but no rule is known!");
						++m_text_iter;
						return process_token();
					}
					// process alphas and digits
					while (!std::isspace(*m_text_iter) && !valid_punctuation()
						&& m_text_iter != m_text_end)
					{
						if (std::isalpha(*m_text_iter))
						{
							internal_token = TOKEN_TYPE::ALPHA;
						}
						else if (std::isdigit(*m_text_iter))
						{
							internal_token = TOKEN_TYPE::NUMERIC;
						}
						else
							throw ParserException("Invalid character.");
						++m_text_iter;
					}

					return process_token();
				}
				/** Getter for the full string being parsed */
				std::string get_full_string() const noexcept
				{
					return std::string(m_text_beginning, m_text_end);
				}

			private:
				/** Iterates one position over the string */
				void advance() { std::next(m_text_end, 1); }
				/**
				 * Processes a given Token. Figures out
				 * what the whole value is and what the type
				 * is.
				 *
				 * @return the current Token
				 */
				Token process_token()
				{
					std::string value;
					switch (internal_token)
					{
					case TOKEN_TYPE::NONE:
						throw ParserException("Unknown character.");
					case TOKEN_TYPE::ALPHA:
					{
						std::copy(m_text_start, m_text_iter, std::back_inserter(value));
						if (str_tolower(value) == "and")
							internal_token = TOKEN_TYPE::AND;
						if (str_tolower(value) == "or")
							internal_token = TOKEN_TYPE::OR;
						if (str_tolower(value) == "to")
							internal_token = TOKEN_TYPE::RANGE;
					}
					break;
					case TOKEN_TYPE::NUMERIC:
					{
						std::copy(m_text_start, m_text_iter, std::back_inserter(value));
						if (!is_valid_numeric(value))
						{
							// if it is not all numeric it is parsed as alpha
							internal_token = TOKEN_TYPE::ALPHA;
						}
					}
					break;
					case TOKEN_TYPE::RPAREN:
					case TOKEN_TYPE::LPAREN:
					{
						std::copy(m_text_start, m_text_iter, std::back_inserter(value));
					}
					break;
					default:
						throw ParserException("tokenizer error");
					}
					return Token(
						internal_token, value, std::distance(m_text_beginning, m_text_start));
				}
				/** Internal function to check if a string represent a digit */
				bool is_valid_numeric(const std::string& token) const noexcept
				{
					for (const char& character : token)
					{
						if (!std::isdigit(character))
							return false;
					}
					return true;
				}
				/** Internal function to check if the punctuation is valid in the DSL */
				bool valid_punctuation() const noexcept
				{
					if (std::ispunct(*m_text_iter))
						if (*m_text_iter == '(' or *m_text_iter == ')')
							return true;
					return false;
				}
				/** The iterator recorded at the start of the string */
				std::string::const_iterator m_text_beginning;
				/** The iterator recorded at the start of the current Token */
				std::string::const_iterator m_text_start;
				/** The current iterator */
				std::string::const_iterator m_text_iter;
				/** The end position of the string */
				std::string::const_iterator m_text_end;
				/** Stores the token type */
				TOKEN_TYPE internal_token;
			};

			/**
			 * The Parser builds the AST with the Tokens
			 * from the Lexer
			 */
			class Parser
			{
			public:
				/**
				 * The parser constructor with a lexer shared pointer
				 * @param lexer the lexer that provides the tokens.
				 */
				Parser(const std::shared_ptr<Lexer>& lexer)
					: m_lexer(lexer)
				{
					try
					{
						m_current_token = m_lexer->get_next_token();
					}
					catch (const std::exception& exc)
					{
						throw exc;
					}
				}

				/**
				 * Parses the sentence being tokenised by the lexer.
				 *
				 * @return the top Node of the AST
				 */
				std::shared_ptr<Node> parse() { return expr(); }

			private:
				/**
				 * Checks TOKEN_TYPE and then iterates.
				 * @param token_type
				 */
				void iter(TOKEN_TYPE token_type)
				{
					if (m_current_token.get_type() == token_type)
					{
						m_current_token = m_lexer->get_next_token();
					}
					else
						throw ParserException(
							"Invalid syntax. Expected type {} but got {} (value: {})",
							get_string_from_enum(token_type),
							get_string_from_enum(m_current_token.get_type()),
							m_current_token.get_value());
				}
				/**
				 * Traceback the error and raise error with a formatted
				 * message, where additional arguments can be passed here.
				 * It will use the current token to traceback.
				 *
				 * @tparam Args types of argument pack
				 * @param error_message the main error message
				 * @param args additional values required in the error message
				 */
				template <typename... Args>
				void traceback(const std::string& error_message, Args... args) const
				{
					traceback(m_current_token, error_message, args...);
				}

				/**
				 * Traceback the error with given token and raise error message.
				 *
				 * @param token that raised the error
				 * @param error_message the main error message
				 */
				void traceback(const Token& token, const std::string& error_message) const
				{
					traceback(token, error_message, "");
				}
				/**
				 * Traceback the error with given token and raise custom error
				 * message, where additional arguments can be passed here.
				 *
				 * @tparam Args types of argument pack
				 * @param token that raised the error
				 * @param error_message the main error message
				 * @param args additional values required in the error message
				 */
				template <typename... Args>
				void traceback(
					const Token& token, const std::string& error_message, Args... args) const
				{
					auto error_msg_offset = ParserException::m_error_msg_start.size();
					auto expression = std::string(error_msg_offset, ' ') + "\""
						+ m_lexer->get_full_string() + "\"";
					std::string error_location(expression.size() + error_msg_offset, ' ');
					auto error_start = token.get_position() + error_msg_offset;
					auto error_end = error_start + token.get_value().size() + 1;

					for (int i = 0; i < expression.size() + error_msg_offset; ++i)
					{
						if (i > error_start && i < error_end)
							error_location[i] = '^';
					}

					expression.append("\n");
					std::string new_error_message
						= error_message + "\n" + expression + error_location;
					throw ParserException(new_error_message, args...);
				}

				/**
				 * Evaluates everythin between logical
				 * operators and decides how to interpret
				 * keywords and whether they should be words
				 * or numbers
				 *
				 * @return the node
				 */
				std::shared_ptr<Node> factor()
				{
					if (m_current_token.get_type() == TOKEN_TYPE::ALPHA)
					{
						TOKEN_TYPE node_type;

						if (str_tolower(m_current_token.get_value()) == "atom")
						{
							node_type = TOKEN_TYPE::ATOM_NAME;
						}
						else if (str_tolower(m_current_token.get_value()) == "residue")
						{
							node_type = TOKEN_TYPE::RESIDUE_NAME;
						}
						else if (str_tolower(m_current_token.get_value()) == "chain")
						{
							node_type = TOKEN_TYPE::CHAIN_NAME;
						}
						else
						{
							traceback("Unknown keyword '{}'.", m_current_token.get_value());
						}
						iter(TOKEN_TYPE::ALPHA);
						if (m_current_token.get_type() == TOKEN_TYPE::NUMERIC)
						{
							// if it is numeric we get the next enum value, which
							// switches alpha to numeric
							++node_type;
						}
						auto node = std::make_shared<Node>(Token(node_type,
							m_current_token.get_value(), m_current_token.get_position()));
						if (m_current_token.get_type() == TOKEN_TYPE::ALPHA)
							iter(TOKEN_TYPE::ALPHA);
						else if (m_current_token.get_type() == TOKEN_TYPE::NUMERIC)
							iter(TOKEN_TYPE::NUMERIC);
						return node;
					}
					else if (m_current_token.get_type() == TOKEN_TYPE::LPAREN)
					{
						auto copy_paren = Token(m_current_token);
						iter(TOKEN_TYPE::LPAREN);
						auto node = expr();
						std::cout << "Solving PAREN: " << node->get_repr() << "\n";
						try
						{
							iter(TOKEN_TYPE::RPAREN);
						}
						catch (const ParserException& e)
						{
							traceback(copy_paren, "Unbalanced paranthesis.");
						}
						return node;
					}
					else
					{
						throw ParserException("Unknown token: {}", m_current_token.get_repr());
					}
				}

				/**
				 * Evaluates an expression.
				 *
				 * @return the top Node
				 */
				std::shared_ptr<Node> expr()
				{
					auto node = factor();
					while (m_current_token.get_type() == TOKEN_TYPE::OR
						|| m_current_token.get_type() == TOKEN_TYPE::AND)
					{
						auto op = Token(m_current_token);
						if (m_current_token.get_type() == TOKEN_TYPE::OR)
							iter(TOKEN_TYPE::OR);
						else if (m_current_token.get_type() == TOKEN_TYPE::AND)
							iter(TOKEN_TYPE::AND);
						auto right = factor();
						node = std::make_shared<Node>(node, op, right);
					}
					if (m_current_token.get_type() != TOKEN_TYPE::EOF_)
						traceback("Unexpected keyword: '{}'.", m_current_token.get_value());
					return node;
				}

				/** The lexer object that returns the tokens */
				std::shared_ptr<Lexer> m_lexer;
				/** The current token returned by the lexer */
				Token m_current_token;
			};
		}

		/**
		 * The DSL (domain specific language) interpreter
		 * of ProStruct.
		 * Rules:
		 *  * 'and' and 'or' are the logical operators
		 *  * 'atom', 'residue' and 'chain' are the keywords
		 *  * '(' and ')' determine resolution order
		 *  * keywords can be either numeric or string
		 *  type which is deduced by the parser.
		 *     * a name is deduced as string if it has any alpha
		 *  	 characters [a-z|A-Z]
		 *     * a name is deduced as numeric if it only has digits
		 *  * whitespaces are ignored and are only used for delimiting
		 *  * keywords are case insensitive
		 *
		 * 	@example
		 * 	- DSLInterpreter("atom CA") returns all CA atoms
		 * 	- DSLInterpreter("residue ALA") returns all alanine residues
		 * 	- DSLInterpreter("atom CA and residue ALA") returns all CA atoms of alanine
		 * residues
		 *
		 */
		class DSLInterpreter
		{
		public:
			DSLInterpreter(const std::string& expression)
			{
				auto lexer = std::make_shared<detail::Lexer>(expression);
				m_parser = std::make_shared<detail::Parser>(lexer);
				// do parsing here so can potentially reuse node
				m_node = m_parser->parse();
			}
			/**
			 * Interpretes the expression and returns
			 * the atom selection.
			 *
			 * @tparam T the float precision type of StructBase
			 * @param m_structure the structure to evaluate
			 * @return the indices selected of the selected atoms
			 */
			template <typename T>
			arma::Col<arma::uword> interpret(const StructBase<T>& structure) const
			{
				arma::Col<arma::uword> max_result(structure.n_atoms());
				arma::uword idx = 0;
				arma::uword i = 0;
				for (const auto& atom : structure.get_atoms())
				{
					if (visit_op(m_node, atom, idx))
					{
						max_result(i) = idx;
						++i;
					}
					++idx;
				}
				return max_result.head(i);
			}

		private:
			/**
			 * Does the atom wise AST evaluation
			 *
			 * @tparam T the Atom type (float or double)
			 * @param node the top Node of the AST
			 * @param atom the atom object to evaluate
			 * @param pos the atom position
			 * @return result of the AST evaluation
			 */
			template <typename T>
			bool visit_op(const std::shared_ptr<detail::Node>& node,
				const std::shared_ptr<Atom<T>>& atom, arma::uword pos) const
			{
				if (node->get_right() == nullptr && node->get_left() == nullptr)
				{
					switch (node->get_type())
					{
					case detail::TOKEN_TYPE::ATOM_NAME:
						return atom->get_name() == node->get_value();
					case detail::TOKEN_TYPE::ATOM_NUMBER:
						return pos == std::stoi(node->get_value());
					case detail::TOKEN_TYPE::RESIDUE_NAME:
						return atom->get_residue_name() == node->get_value();
					case detail::TOKEN_TYPE::RESIDUE_NUMBER:
						return atom->get_residue_number() == node->get_value();
					case detail::TOKEN_TYPE::CHAIN_NAME:
						return atom->get_residue().get_chain_name() == node->get_value();
					default:
						throw detail::ParserException("Unknown keyword: {}", node->get_repr());
					}
				}
				switch (node->get_op_type())
				{
				case detail::TOKEN_TYPE::AND:
					return visit_op(node->get_left(), atom, pos)
						&& visit_op(node->get_right(), atom, pos);
				case detail::TOKEN_TYPE::OR:
					return visit_op(node->get_left(), atom, pos)
						|| visit_op(node->get_right(), atom, pos);
				default:
					throw detail::ParserException(
						"Unknown op: {}", detail::get_string_from_enum(node->get_op_type()));
				}
			}
			/** The parser instance for given expression */
			std::shared_ptr<detail::Parser> m_parser;
			/** The top node of the AST */
			std::shared_ptr<detail::Node> m_node;
		};
	}
}

#endif // PROSTRUCT_AST_H
