solution "default"
  configurations { "Debug", "Release" }

  project "lumi_analysis"
    language "C++"
    kind     "ConsoleApp"
    files  { "include/*.h", "src/*.cc", "main/*.cc" }
    includedirs { "include" }

    configuration { "Debug*" }
      buildoptions { "-std=c++0x", "$(shell root-config --cflags)" }
      linkoptions { "$(shell root-config --libs)", "-lMinuit" }
      defines { "_DEBUG", "DEBUG" }
      flags   { "Symbols", "ExtraWarnings" }

    configuration { "Release*" }
      buildoptions { "-std=c++0x", "$(shell root-config --cflags)" }
      linkoptions { "$(shell root-config --libs)", "-lMinuit" }
      defines { "NDEBUG" }
      flags   { "Optimize" }
