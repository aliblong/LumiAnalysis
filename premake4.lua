local priv = "${HOME}/privatemodules/"

solution "default"
  configurations { "Debug", "Release" }

  project "lumi_analysis"
    language "C++"
    kind     "ConsoleApp"
    files  { "include/*.h", "src/*.cc", "main/*.cc" }
    includedirs { "include", priv.."boost_1_57_0", priv.."expected/include", priv.."Mach7/code", priv.."root/include" }

    configuration { "Debug*" }
      buildoptions { "-std=c++1y", "-Wno-unused-local-typedefs", "-pthread", "-m64" } -- "-Wno-deprecated-declarations", "$(shell root-config --cflags)" }--, "-fopenmp" }
      linkoptions { "$(shell root-config --libs)", "-lMinuit" } --, priv.."Mach7/code" }--, "-fopenmp" }
      defines { "_DEBUG", "DEBUG" }
      flags   { "Symbols", "ExtraWarnings" }

    configuration { "Release*" }
      buildoptions { "-std=c++1y", "$(shell root-config --cflags)" }
      linkoptions { "$(shell root-config --libs)", "-lMinuit" }
      defines { "NDEBUG" }
      flags   { "Optimize" }
