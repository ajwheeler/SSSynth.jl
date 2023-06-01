@testset "fit EWs" begin
    atm = read_model_atmosphere("data/sun.mod")

    Feline = Korg.Line(5e-5, 0.0, Korg.species"Fe I", 0.0)
    FeHline = Korg.Line(5.5e-5, 0.0, Korg.species"FeH", 0.0)
    Oline = Korg.Line(7e-5, 0.0, Korg.species"O I", 0.0)

    A_X = format_A_X()

    @test_throws ArgumentError Korg.fit_EWs(atm, [Feline, Feline], A_X, [0,0])
    @test_throws ArgumentError Korg.fit_EWs(atm, [FeHline], A_X, [0,0])
    @test_throws ArgumentError Korg.fit_EWs(atm, [Feline, Oline], A_X, [0,0])
end