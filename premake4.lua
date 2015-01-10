local priv = "${HOME}/privatemodules/"

solution "default"
  configurations { "Debug", "Release" }

  project "lumi_analysis"
    language "C++"
    kind     "ConsoleApp"
    files  { "include/*.h", "src/*.cc", "main/*.cc" }
    includedirs { "include", priv.."boost_1_57_0", priv.."expected/include", priv.."Mach7/code" }

    configuration { "Debug*" }
      buildoptions { "-std=c++11", "-Wno-unused-local-typedefs", "$(shell root-config --cflags)" }--, "-fopenmp" }
      linkoptions { "$(shell root-config --libs)", "-lMinuit" } --, priv.."Mach7/code" }--, "-fopenmp" }
      defines { "_DEBUG", "DEBUG" }
      flags   { "Symbols", "ExtraWarnings" }

    configuration { "Release*" }
      buildoptions { "-std=c++11", "$(shell root-config --cflags)" }
      linkoptions { "$(shell root-config --libs)", "-lMinuit" }
      defines { "NDEBUG" }
      flags   { "Optimize" }
