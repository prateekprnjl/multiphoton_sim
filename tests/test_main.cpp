#define CATCH_CONFIG_RUNNER

#include <catch2/catch_session.hpp>
#include <petscsys.h>

int main(int argc, char* argv[])
{
    PetscInitialize(&argc, &argv, nullptr, "Running tests.");

    int result = Catch::Session().run(argc, argv);

    PetscFinalize();

    return result;
}
