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
#include <iostream>
#include <unordered_map>
#include <variant>

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
			class ParserException : public std::exception
			{
			public:
				ParserException(const std::string& message)
					: m_message(message)
				{
				}

				const char* what() const noexcept
				{
					return format(fmt("Parser error: {}"), m_message).c_str();
				}

			private:
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
			inline std::string get_string_from_enum(TOKEN_TYPE token_type)
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
			}

			/**
			 * Convert a string to lower case.
			 * Straight from https://en.cppreference.com/w/cpp/string/byte/tolower.
			 * @param s
			 * @return
			 */
			inline std::string str_tolower(std::string s)
			{
				std::transform(s.begin(), s.end(), s.begin(),
					[](unsigned char c) { return std::tolower(c); }
				);
				return s;
			}

			/**
			 * A token is a combination of a key that helps
			 * identify the type of the value of the Token.
			 */
			class Token
			{
			public:
				Token()
					: m_type(TOKEN_TYPE::NONE)
					, m_value("")
				{
				}

				Token(TOKEN_TYPE name, const std::string& value)
					: m_type(name)
					, m_value(value)
				{
				}

				TOKEN_TYPE get_type() const noexcept { return m_type; }

				std::string get_value() const noexcept { return m_value; }

				void print() const noexcept { std::cout << get_repr() << "\n"; }

				std::string get_repr() const noexcept
				{
					return format(
						fmt("Token(value={}, type={})"), m_value, get_string_from_enum(m_type));
				}

			private:
				TOKEN_TYPE m_type;
				std::string m_value;
			};

			/**
			 * The Node in the AST tree with various
			 * helper functions.
			 */
			class Node
			{
			public:
				Node() {}

				Node(Token token)
					: m_token(token)
				{
				}

				Node(const std::shared_ptr<Node>& left, Token token,
					const std::shared_ptr<Node>& right)
					: m_left(left)
					, m_token(token)
					, m_right(right)
				{
					if (m_token.get_type() != TOKEN_TYPE::OR
						&& m_token.get_type() != TOKEN_TYPE::AND)
						throw ParserException(
							format(fmt("Expected op to be either of type {} or {}, but got {}"),
								get_string_from_enum(TOKEN_TYPE::OR),
								get_string_from_enum(TOKEN_TYPE::AND),
								get_string_from_enum(m_token.get_type())));
				}

				TOKEN_TYPE get_op_type() const noexcept { return m_token.get_type(); }

				std::string get_right_repr() const noexcept
				{
					if (m_right)
						return get_right()->get_repr();
					else
						return std::string("N/A");
				}

				std::string get_left_repr() const noexcept
				{
					if (m_right)
						return get_left()->get_repr();
					else
						return std::string("N/A");
				}

				std::string get_repr() const noexcept
				{
					return format(fmt("Node(left={}, right={}, token={})"), get_left_repr(),
						get_right_repr(), m_token.get_repr());
				}

				std::shared_ptr<Node> get_left() const noexcept { return m_left; }

				std::shared_ptr<Node> get_right() const noexcept { return m_right; }

				TOKEN_TYPE get_type() const noexcept { return m_token.get_type(); }

				std::string get_value() const noexcept { return m_token.get_value(); }

			private:
				std::shared_ptr<Node> m_left;
				std::shared_ptr<Node> m_right;
				Token m_token;
			};

			/**
			 * The Lexer parses a string and represents with
			 * Tokens. The Lexer ignores whitespaces.
			 */
			class Lexer
			{
			public:
				Lexer(const std::string& text)
					: m_text_start(text.cbegin())
					, m_text_iter(text.cbegin())
					, m_text_end(text.cend())
					, internal_token(TOKEN_TYPE::NONE)
					, m_freeze(false)
				{
				}

				void advance() { std::next(m_text_end, 1); }

				Token process_token()
				{
					std::string value;
					switch (internal_token)
					{
					case TOKEN_TYPE::NONE:
					{
						++m_text_iter;
						return get_next_token();
					}
					break;
					case TOKEN_TYPE::ALPHA:
					{
						std::copy(m_text_start, m_text_iter, std::back_inserter(value));
						if (value == "and")
							internal_token = TOKEN_TYPE::AND;
						if (value == "or")
							internal_token = TOKEN_TYPE::OR;
					}
					break;
					case TOKEN_TYPE::NUMERIC:
					{
						std::copy(m_text_start, m_text_iter, std::back_inserter(value));
					}
					break;
					case TOKEN_TYPE::RPAREN:
					case TOKEN_TYPE::LPAREN:
					{
						value = std::string(m_text_iter, m_text_iter + 1);
					}
					break;
					default:
						throw ParserException("tokenizer error");
					}
					if (!m_freeze)
						++m_text_iter;
					return Token(internal_token, value);
				}

				Token get_next_token()
				{
					bool look_ahead = false;

					while (m_text_iter != m_text_end)
					{
						if (std::isspace(*m_text_iter))
						{
							look_ahead = false;
							m_freeze = false;
							auto result = process_token();
							internal_token = TOKEN_TYPE::NONE;
							return result;
						}
						else if (std::isalpha(*m_text_iter))
						{
							if (!look_ahead)
								m_text_start = m_text_iter;
							look_ahead = true;
							m_freeze = false;
							internal_token = TOKEN_TYPE::ALPHA;
							++m_text_iter;
						}
						else if (std::isdigit(*m_text_iter))
						{
							if (!look_ahead)
								m_text_start = m_text_iter;
							// if the token started as a alpha it will stay an alpha
							if (internal_token == TOKEN_TYPE::ALPHA)
								internal_token = TOKEN_TYPE::ALPHA;
							else
								internal_token = TOKEN_TYPE::NUMERIC;
							look_ahead = true;
							m_freeze = false;
							++m_text_iter;
						}
						else if (std::ispunct(*m_text_iter))
						{
							look_ahead = false;
							if (*m_text_iter == '(')
							{
								m_freeze = false;
								internal_token = TOKEN_TYPE::LPAREN;
							}
							else if (*m_text_iter == ')')
							{
								if (m_freeze)
									m_freeze = false;
								else
								{
									m_freeze = true;
									return process_token();
								}
								internal_token = TOKEN_TYPE::RPAREN;
							}
							else
								throw ParserException(
									format(fmt("Unknown punctuation character {}"), *m_text_iter));
							auto result = process_token();
							internal_token = TOKEN_TYPE::NONE;
							return result;
						}
						else
							throw ParserException(
								format(fmt("Unknown character {}"), *m_text_iter));
					}
					if (!m_freeze)
					{
						m_freeze = true;
						return process_token();
					}
					else
						return Token(TOKEN_TYPE::EOF_, "EOF");
				}

			private:
				std::string::const_iterator m_text_start;
				std::string::const_iterator m_text_iter;
				std::string::const_iterator m_text_end;
				TOKEN_TYPE internal_token;
				bool m_freeze;
			};

			/**
			 * The Parser builds the AST with the Tokens
			 * from the Lexer
			 */
			class Parser
			{
			public:
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

				std::shared_ptr<Node> parse() {return expr();}
			private:
				/**
				 * Checks TOKEN_TYPE and then iterates.
				 * @param token_type
				 */
				void iter(TOKEN_TYPE token_type)
				{
					// std::cout << "Iterating from " << m_current_token.get_repr() << "\n";
					if (m_current_token.get_type() == TOKEN_TYPE::EOF_)
						return;
					else if (m_current_token.get_type() == token_type)
					{
						m_current_token = m_lexer->get_next_token();
						// std::cout << "New token: " << m_current_token.get_repr() << "\n";
					}
					else
						throw ParserException(
							format(fmt("Invalid syntax. Expected type {} but got {} (value: {})"),
								get_string_from_enum(token_type),
								get_string_from_enum(m_current_token.get_type()),
								m_current_token.get_value()));
				}

				std::shared_ptr<Node> factor()
				{
					// std::cout << "factoring: " << m_current_token.get_repr() << "\n";
					if (m_current_token.get_type() == TOKEN_TYPE::ALPHA)
					{
						// std::cout << "it's an alpha\n";
						// std::cout << "factoring: " << m_current_token.get_repr() << "\n";
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
							node_type = m_current_token.get_type();
						}
						iter(TOKEN_TYPE::ALPHA);
						if (m_current_token.get_type() == TOKEN_TYPE::NUMERIC)
						{
							++node_type; // if it is numeric we get the next enum value, which
										 // switches alpha to numeric
						}
						auto node
							= std::make_shared<Node>(Token(node_type, m_current_token.get_value()));
						// std::cout << "factored: " << node->get_repr() << "\n";
						if (m_current_token.get_type() == TOKEN_TYPE::ALPHA)
							iter(TOKEN_TYPE::ALPHA);
						else if (m_current_token.get_type() == TOKEN_TYPE::NUMERIC)
							iter(TOKEN_TYPE::NUMERIC);
						return node;
					}
					else if (m_current_token.get_type() == TOKEN_TYPE::LPAREN)
					{
						iter(TOKEN_TYPE::LPAREN);
						auto node = term();
						// std::cout << node->get_repr() << "\n";
						iter(TOKEN_TYPE::RPAREN);
						return node;
					}
					else
					{
						throw ParserException(
							format(fmt("Unknown token: {}"), m_current_token.get_repr()));
						// std::cout << format(fmt("Unknown token: {}"), m_current_token.get_repr())
						// << "\n"; iter(m_current_token.get_type()); std::cout << "Skipped\n";
						// std::cout << format(fmt("Next token: {}"), m_current_token.get_repr()) <<
						// "\n"; return nullptr;
					}
				}

				std::shared_ptr<Node> term()
				{
					auto node = factor();
					if (m_current_token.get_type() == TOKEN_TYPE::ALPHA)
					{
						// std::cout << "op: " << std::make_shared<Token>(m_current_token)->get_repr() << "\n";
						auto op = Token(m_current_token);
						iter(m_current_token.get_type());
						auto right = factor();
						//			std::cout << "right: " << right->get_repr() << "\n";
						//			std::cout << "making node\n";
						return std::make_shared<Node>(node, op, right);
					}
					else
						return node;
				}

				/**
				 * Evaluates an expression.
				 *
				 */
				std::shared_ptr<Node> expr()
				{
					std::shared_ptr<Node> node = term();
					while (m_current_token.get_type() == TOKEN_TYPE::OR
						|| m_current_token.get_type() == TOKEN_TYPE::AND)
					{
						// std::cout << "Current token: " << m_current_token.get_repr() << "\n";
						// std::cout << "CURRENT LEFT NODE: " << node->get_repr() << "\n";
						auto maybe_op_token = Token(m_current_token);
						// std::cout << "CURRENT TOKEN: " << maybe_op_token.get_repr() << "\n";
						if (m_current_token.get_type() == TOKEN_TYPE::OR)
							iter(TOKEN_TYPE::OR);
						if (m_current_token.get_type() == TOKEN_TYPE::AND)
							iter(TOKEN_TYPE::AND);
						auto right = term();
						// if (right)
						// {
						// 	std::cout << "CURRENT RIGHT NODE: " << right->get_repr() << "\n";
						// }
						node = std::make_shared<Node>(node, maybe_op_token, right);
						
						// std::cout << "End node: " << node->get_repr() << "\n";
					}
					return node;
				}
				std::shared_ptr<Lexer> m_lexer;
				Token m_current_token;
			};
		}

		/**
		 * The DSL (domain specific language) interpreter
		 * of ProStruct.
		 * Rules:
		 * 	- 'and' and 'or' are the logical operators
		 * 	- 'atom', 'residue' and 'chain' are the keywords
		 * 	- '(' and ')' determine resolution order
		 * 	- keywords can either be either of numeric or string
		 * 	type which is deduced by the parser.
		 * 	- whitespaces are ignored and are only used for delimiting
		 * 	- keywords are case insensitive
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
			}

			template <typename T>
			arma::Col<arma::uword> interpret(const StructBase<T>& m_structure) const
			{
				arma::Col<arma::uword> max_result(m_structure.n_atoms());
				auto node = m_parser->parse();
				arma::uword idx = 0;
				arma::uword i = 0;
				for (const auto& atom : m_structure.get_atoms())
				{
					if (visit_op(node, atom, idx))
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
				const std::shared_ptr<Atom<T>>& atom, int pos) const
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
						throw detail::ParserException(
							format(fmt("Unknown keyword: {}"), node->get_repr()));
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
					throw detail::ParserException(format(
						fmt("Unknown op: {}"), detail::get_string_from_enum(node->get_op_type())));
				}
			}

			std::shared_ptr<detail::Parser> m_parser;
		};
	}
}

#endif // PROSTRUCT_AST_H
