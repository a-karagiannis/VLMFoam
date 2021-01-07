Info<< "Reading NucleationCoefficients\n" << endl;

IOdictionary NucleationCoefficients
(
    IOobject
    (
        "NucleationCoefficients",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const scalar alpha_coeff(readScalar(NucleationCoefficients.subDict("thermalCorrections").lookup("alpha_coeff")));
const scalar q_c(readScalar(NucleationCoefficients.subDict("thermalCorrections").lookup("q_c")));

