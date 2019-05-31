#include <prostruct/prostruct.h>
#include <prostruct/utils/ast.h>

using namespace prostruct;
using namespace prostruct::parser::detail;

int main()
{
	auto test = std::string("residue 10 to 20 and atom CA");
	auto pdb = PDB<float>("./build/test.pdb");
	// std::cout << pdb.n_atoms() << "\n";
	// auto test = std::string("atom CA or atom CB and residue 10");
	auto test_parser = std::string(test);
	auto lexer = std::make_shared<Lexer>(test_parser);
	int i = 0;
	while (i<10)
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
		{
			std::cout << "END\n";
			break;
		}
		i++;
	}
	std::cout << "\n=====\n";
	std::shared_ptr<Node> node;
	lexer = std::make_shared<Lexer>(test);
	auto parser = Parser(lexer);
	try
	{
		node = parser.parse();
	}
	catch (const std::exception& exc)
	{
		std::cout << exc.what() << "\n" << "ABORTING\n";
		return 0;
	}
	std::cout << "RESULT: " << node->get_repr() << "\n";
	auto dsl_interp = prostruct::parser::DSLInterpreter(test);
	dsl_interp.interpret(pdb).print();
}
