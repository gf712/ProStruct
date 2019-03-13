#include <prostruct/utils/ast.h>
#include <prostruct/prostruct.h>

using namespace prostruct::parser::detail;

int main()
{
	// auto test = std::string("(atom CA or atom CB) and residue 10");
	auto pdb = PDB<float>("./build/test.pdb");
	std::cout << pdb.n_atoms() << "\n";
	// auto test = std::string("atom CA or atom CB and residue 10");
	auto test = std::string("atom CA");
	auto lexer = std::make_shared<Lexer>(test);
	while (true)
	{
		Token token;
		try {
			token = lexer->get_next_token();
		}
		catch (const std::exception& exc)
		{
			std::cout << exc.what() << "\n" << "ABORTING\n";
			return 0;
		}
		token.print();
		if (token.get_type() == TOKEN_TYPE::EOF_)
			break;
	}

	// auto node = Node(std::make_shared<Keyword>(Token(TOKEN_TYPE::ALPHA, "atom"), Token(TOKEN_TYPE::ALPHA, "CA")));
	// std::cout << node.get_repr() << "\n";
	lexer = std::make_shared<Lexer>(test);
	auto parser = Parser(lexer);
	auto node = parser.expr();
	std::cout << node->get_repr() << "\n";
	auto dsl_interp = prostruct::parser::DSLInterpreter("atom CA");
	dsl_interp.interpret(pdb).print();
}
						
	// 					BinaryOp
	// 						|				
	// 			|						|						
	// 		BinaryOp				residue 10
	// 			|		
	// 	|				|
	// atom CA 		  atom CB