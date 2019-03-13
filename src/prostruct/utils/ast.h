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

			enum class TOKEN_TYPE
			{
				NONE,
				NUMERIC,
				ALPHA,
				RPAREN,
				LPAREN,
				AND,
				OR,
				EOF_
			};

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
				}
			}

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
						fmt("Token(value={}, type={})"), m_value, static_cast<int>(m_type));
				}

			private:
				TOKEN_TYPE m_type;
				std::string m_value;
			};

			class Keyword
			{
			public:
				Keyword() {}

				Keyword(const Token& type, const Token& value)
					: m_type(type)
					, m_value(value)
				{
					if (type.get_type() != TOKEN_TYPE::ALPHA)
						throw ParserException(
							format(fmt("A keyword type must be an alpha, but got {} (value: {})"),
								static_cast<int>(type.get_type()), type.get_value()));
					if (type.get_type() != TOKEN_TYPE::ALPHA
						&& type.get_type() != TOKEN_TYPE::NUMERIC)
						throw ParserException(format(fmt("A keyword value must be an alpha or "
														 "numeric, but got {} (value: {})"),
							static_cast<int>(type.get_type()), type.get_value()));
				}

				void print() const noexcept { std::cout << get_repr() << "\n"; }

				std::string get_repr() const noexcept
				{
					return format(
						fmt("Keyword(value={}, type={})"), m_value.get_value(), m_type.get_value());
				}

				Token get_type() const noexcept { return m_type; }

				Token get_value() const noexcept { return m_value; }

			private:
				Token m_type;
				Token m_value;
			};

			class Node
			{
			public:
				Node() {}

				Node(const std::variant<std::shared_ptr<Keyword>, std::shared_ptr<Node>>& left)
					: m_left(left)
				{
				}

				Node(const std::variant<std::shared_ptr<Keyword>, std::shared_ptr<Node>>& left,
					const std::shared_ptr<Token>& op,
					const std::variant<std::shared_ptr<Keyword>, std::shared_ptr<Node>>& right)
					: m_left(left)
					, m_op(op)
					, m_right(right)
				{
					if (op->get_type() != TOKEN_TYPE::OR && op->get_type() != TOKEN_TYPE::AND)
						throw ParserException(
							format(fmt("Expected op to be either of type {} or {}, but got {}"),
								static_cast<int>(TOKEN_TYPE::OR), static_cast<int>(TOKEN_TYPE::AND),
								static_cast<int>(m_op->get_type())));
				}

				TOKEN_TYPE get_op_type() const noexcept { return m_op->get_type(); }

				std::string get_right_repr() const noexcept
				{
					return std::visit(
						[](auto&& arg) {
							if (arg != nullptr)
								return arg->get_repr();
							else
								return std::string("N/A");
						},
						m_right);
				}

				template <typename T>
				bool eval_left(const std::shared_ptr<Atom<T>>)
				{
					TOKEN_TYPE token_type
						= std::visit([](auto&& arg) { return arg->get_type(); }, get_left());

					switch (token_type)
					{
					case TOKEN_TYPE::ALPHA:
					{
					}
					}
				}

				template <typename T>
				bool eval_right(const std::shared_ptr<Atom<T>>)
				{
				}

				std::string get_left_repr() const noexcept
				{
					return std::visit(
						[](auto&& arg) {
							if (arg != nullptr)
								return arg->get_repr();
							else
								return std::string("N/A");
						},
						m_left);
				}

				std::string get_op_repr() const noexcept
				{
					if (m_op != nullptr)
						return m_op->get_repr();
					else
						return std::string("N/A");
				}

				std::string get_repr() const noexcept
				{
					return format(fmt("Node(left={}, op={}, right={})"), get_left_repr(),
						get_op_repr(), get_right_repr());
				}

				std::variant<std::shared_ptr<Keyword>, std::shared_ptr<Node>> get_left() const
					noexcept
				{
					return m_left;
				}

				std::variant<std::shared_ptr<Keyword>, std::shared_ptr<Node>> get_right() const
					noexcept
				{
					return m_right;
				}

			private:
				std::variant<std::shared_ptr<Keyword>, std::shared_ptr<Node>> m_left;
				std::shared_ptr<Token> m_op;
				std::variant<std::shared_ptr<Keyword>, std::shared_ptr<Node>> m_right;
			};

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
						value = "";
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
						return Token(TOKEN_TYPE::EOF_, "");
				}

			private:
				std::string::const_iterator m_text_start;
				std::string::const_iterator m_text_iter;
				std::string::const_iterator m_text_end;
				TOKEN_TYPE internal_token;
				bool m_freeze;
			};

			/**
			 * Parser to build an AST to find atoms with specific
			 * properties.
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

				/**
				 * Checks TOKEN_TYPE and then iterates.
				 * @param token_type
				 */
				void iter(TOKEN_TYPE token_type)
				{
					if (m_current_token.get_type() == TOKEN_TYPE::NONE)
					{
						m_current_token = m_lexer->get_next_token();
						iter(token_type);
					}
					else if (m_current_token.get_type() == token_type)
						m_current_token = m_lexer->get_next_token();
					else
						throw ParserException(
							format(fmt("Invalid syntax. Expected type {} but got {} (value: {})"),
								get_string_from_enum(token_type),
								get_string_from_enum(m_current_token.get_type()),
								m_current_token.get_value()));
				}

				std::shared_ptr<Node> factor()
				{
					//				std::cout << "factoring: " << m_current_token.get_repr() <<
					//"\n";
					if (m_current_token.get_type() == TOKEN_TYPE::ALPHA)
					{
						//			std::cout << "it's an alpha\n";
						Token current_token = m_current_token;
						iter(TOKEN_TYPE::ALPHA);
						//			std::cout << "factoring: " << m_current_token.get_repr() <<
						//"\n";
						auto node = std::make_shared<Node>(
							std::make_shared<Keyword>(current_token, m_current_token));
						//					std::cout << "factored: " << node->get_repr() << "\n";
						if (m_current_token.get_type() == TOKEN_TYPE::ALPHA)
							iter(TOKEN_TYPE::ALPHA);
						if (m_current_token.get_type() == TOKEN_TYPE::NUMERIC)
							iter(TOKEN_TYPE::NUMERIC);
						return node;
					}
					else if (m_current_token.get_type() == TOKEN_TYPE::LPAREN)
					{
						iter(TOKEN_TYPE::LPAREN);
						auto node = term();
						//			std::cout << node->get_repr() << "\n";
						iter(TOKEN_TYPE::RPAREN);
						return node;
					}
				}

				std::shared_ptr<Node> term()
				{
					auto node = factor();
					if (m_current_token.get_type() == TOKEN_TYPE::OR
						|| m_current_token.get_type() == TOKEN_TYPE::AND)
					{
						//			std::cout << "it's an OR\n";
						//			std::cout << "left: "
						//					  << std::visit([](auto arg) { return arg->get_repr();
						//}, node->get_left())
						//					  << "\n";
						//			std::cout << "op: " <<
						// std::make_shared<Token>(m_current_token)->get_repr() << "\n";
						auto op = std::make_shared<Token>(m_current_token);
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
					while (m_current_token.get_type() != TOKEN_TYPE::EOF_)
					{
						//					std::cout << "CURRENT LEFT NODE: " << node->get_repr()
						//<<
						//"\n";
						auto maybe_op_token = std::make_shared<Token>(m_current_token);
						if (maybe_op_token->get_type() == TOKEN_TYPE::OR)
							iter(TOKEN_TYPE::OR);
						else if (maybe_op_token->get_type() == TOKEN_TYPE::AND)
							iter(TOKEN_TYPE::AND);
						//			std::cout << "CURRENT TOKEN: " << maybe_op_token->get_repr() <<
						//"\n";
						auto left = term();
						//			std::cout << "CURRENT RIGHT NODE: " << left->get_repr() << "\n";
						node = std::make_shared<Node>(node, maybe_op_token, left);
					}
					return node;
				}

			private:
				std::shared_ptr<Lexer> m_lexer;
				Token m_current_token;
			};
		}

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
				auto node = m_parser->expr();
				std::cout << node->get_repr() << "\n";
				arma::uword idx = 0;
				arma::uword i = 0;
				for (const auto& atom : m_structure.get_atoms())
				{
					if (visit_op(node, atom))
					{
						max_result(i) = idx;
						++i;
					}
					++idx;
				}
				return max_result.head(i);
			}

			template <typename T>
			bool visit_op(
				const std::shared_ptr<detail::Node>& node, const std::shared_ptr<Atom<T>>& atom) const
			{
				if (std::visit([](auto&& arg) { return arg == nullptr; }, node->get_right()))
				{
					auto key = std::get<std::shared_ptr<detail::Keyword>>(node->get_left());
					if (key->get_type().get_value() == "atom")
					{
						return atom->get_name() == key->get_value().get_value();
					}
					else if (key->get_type().get_value() == "residue")
					{
						return atom->get_residue_name() == key->get_value().get_value();
					}
					else
					{
						throw detail::ParserException(
							format(fmt("Unknown keyword: {}"), key->get_repr()));
					}
				}
				switch (node->get_op_type())
				{
				case detail::TOKEN_TYPE::AND:
				{
					return visit_op(std::get<std::shared_ptr<detail::Node>>(node->get_left()), atom)
						&& visit_op(
							std::get<std::shared_ptr<detail::Node>>(node->get_right()), atom);
				}
				case detail::TOKEN_TYPE::OR:
				{
					return visit_op(std::get<std::shared_ptr<detail::Node>>(node->get_left()), atom)
						|| visit_op(
							std::get<std::shared_ptr<detail::Node>>(node->get_right()), atom);
				}
				default:
					throw detail::ParserException(format(
						fmt("Unknown op: {}"), detail::get_string_from_enum(node->get_op_type())));
				}
			}

		private:
			std::shared_ptr<detail::Parser> m_parser;
		};
	}
}

#endif // PROSTRUCT_AST_H
