# The DUNE configuration of the Technical Design Report (DUNE collaboration,
# arXiv:2103.04797). The AEDL files are in examples/DUNE; the data tables are
# the ancillary files of the paper on arXiv, downloaded on first use.

const DUNE_TABLES_URL = "https://arxiv.org/src/2103.04797v1/anc/dune_globes"

# SHA-256 of each table, relative to the directory of DUNE_GLoBES.glb
const DUNE_TABLES = Dict(
    "eff/GLB_app_FHC_eff.dat" => "3120bb60ce5239912474bf336d501bc3d3a6853d49c12283536eeea8c5cbdf7b",
    "eff/GLB_app_RHC_eff.dat" => "71191cff805dfbc3907e491e19cbeb8d77d6e15c651df289564e07b7c3c31408",
    "eff/GLB_dis_FHC_eff.dat" => "6cb7ed079df4db2962f89c0664b7af6a951a00579a395835683a05f4b8c0a18b",
    "eff/GLB_dis_RHC_eff.dat" => "f4f150f15fe8f252ba28be897d7a95e124eaec4431bfd999a079338dde58257a",
    "eff/post_app_FHC_aNC_bkg.txt" => "03a593d195e36f6f2a18fdc777940812799b663942acb630e42155feb54c921b",
    "eff/post_app_FHC_anuNC_bkg.txt" => "55ebbca4b678a265da7a1868f44039a3f1264ab4a87b1bda8998ae41317de6f2",
    "eff/post_app_FHC_NC_bkg.txt" => "7421beabb557dd0931c3473a4f57af9811acb34ea36c55d3fbf57af2313c0f9c",
    "eff/post_app_FHC_nuebar_bkg.txt" => "7bc18248d1072f28b2d20a5952150c421a4311429242779fb733ad3f26cb6b2c",
    "eff/post_app_FHC_nuebar_sig.txt" => "209d6237277a8f5cb6585c8a312cfe2a3613f239c2c065be846fe9fb8cce9b96",
    "eff/post_app_FHC_nue_bkg.txt" => "c9c5dfca5f7020f5dfed44a65c12b2107c33020490d5701846222186a01727ae",
    "eff/post_app_FHC_nue_sig.txt" => "df30f5beb502620acb7d1974f6f5e26aaef9578814e6f59ab86dd31c6d49fbb6",
    "eff/post_app_FHC_numubar_bkg.txt" => "9584acc8bf8656fa57d2d01ecc53b951e545de6a10966c928324620583fb43a3",
    "eff/post_app_FHC_numu_bkg.txt" => "5675aa878b27ee0e192837482f905420962b54ba44a54c0bef2ccb1ef3e9e087",
    "eff/post_app_FHC_nutaubar_bkg.txt" => "3046de33a5dff52e50516d6b0b3fa7848445945069c8ab7eeeb0295f862efd00",
    "eff/post_app_FHC_nutau_bkg.txt" => "d62f3b50afaa8cf1d21903a82eecef2754726ed8c0a232c14476fc37358cc821",
    "eff/post_app_RHC_aNC_bkg.txt" => "67f3e1120d8f8b60768044cd12f74f7ddf6f6771ba3befbe89bc364343a7605c",
    "eff/post_app_RHC_NC_bkg.txt" => "cfd74b0e4088ed8faf50811b7c4ec96dd8279d23f91ad5bea608e26fcdc73227",
    "eff/post_app_RHC_nuebar_bkg.txt" => "c9f1ee516f42d14c524f7b9978055c86d4f2ce9a76132cc0c503a6c25257b041",
    "eff/post_app_RHC_nuebar_sig.txt" => "203c93df7764db59f588a4981efc9376a7fb559a16d2889acb934b04ba2a2725",
    "eff/post_app_RHC_nue_bkg.txt" => "059dcc2c26f3596254bc383516a56fd17dd1e126468ff40039f046d41779a581",
    "eff/post_app_RHC_nue_sig.txt" => "ff2320f4748057484c68af1aa22f4a68bb4dbf879a5323f080fc9202f2f3783f",
    "eff/post_app_RHC_numubar_bkg.txt" => "de295f14c9db34ec8c1b09aa09e7ba8679541b17666d653db812003d52ab80ea",
    "eff/post_app_RHC_numu_bkg.txt" => "83ef2690ba84b4c6d929f508aa0f9377567aa1bc9f10ce9b7802bfbe23aa0f61",
    "eff/post_app_RHC_nutaubar_bkg.txt" => "4074a359863e828d07e85d4e94b48abb0068e1ed7b287046f93751276004a5d7",
    "eff/post_app_RHC_nutau_bkg.txt" => "ae38838c7d9a1552da2bf07a954f1dd643058aa7f4e245a49e6b72b93b2c214d",
    "eff/post_dis_FHC_aNC_bkg.txt" => "3d7e9335aacd3187632c4f2402ff134f257238fad5302127711d9f8704c714a4",
    "eff/post_dis_FHC_NC_bkg.txt" => "4033b070b40142bed92dd01e77e1a4a024a72a563ef97cc87419d3d8a4aaf7a8",
    "eff/post_dis_FHC_numubar_sig.txt" => "a84bce03be38ee6c93e9d301cb8d7879a2eb4ecda86bcbec4e70611b4c1abbb2",
    "eff/post_dis_FHC_numu_sig.txt" => "e3768bf238d87e288b307e54b851110e655f983989c9199c7a1d598fe3691fee",
    "eff/post_dis_FHC_nutaubar_bkg.txt" => "66e8bf79ebd262cd37139889bb693f65783be6e959c78824ae9dc24f4d20637c",
    "eff/post_dis_FHC_nutau_bkg.txt" => "721ed310ede326d1e3aeb3bf01e1d22bde390e6d5cf02886f21605f64210449c",
    "eff/post_dis_RHC_aNC_bkg.txt" => "3b8b1832580963f9495bc13b03956be698b5fcee38a8e85a30482eff6f38b04f",
    "eff/post_dis_RHC_NC_bkg.txt" => "06623fead7931c3c3a61285e342f005185d1110e32c66b2835a00ab8958efc51",
    "eff/post_dis_RHC_numubar_sig.txt" => "4522bf76215f7f3b63ce05749fce901849b24e898fba268b18f7aaa52d4eaa5b",
    "eff/post_dis_RHC_numu_sig.txt" => "0acc732e398c529eed0fc179a1bd384760fcc72ae06fc5068cdc58d23973fcbb",
    "eff/post_dis_RHC_nutaubar_bkg.txt" => "b70e9d668474ef49ace25996434cb1b11f46e6be112bce487a6ea2ede8392ac7",
    "eff/post_dis_RHC_nutau_bkg.txt" => "22569aacf97d4b134a74e0357131ad58718bf2514e11629509a34e05a9e237dd",
    "flux/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_antineutrino_LBNEFD_globes_flux.txt" => "7136b7e9aef4d1c38802032b399aa19b78f1d64c8fe12536c0d3733e6aacd97b",
    "flux/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEFD_globes_flux.txt" => "0dd2eb8132eff8ea7d4bceaa1255a96bc95470262787d0b25e97ff8a21371500",
    "smr/app_aNC_bkg.txt" => "5e3a73a79a183ffe0453c137d1161958bc49acdff7182b38964737a1eb5bd1fa",
    "smr/app_NC_bkg.txt" => "45661d6c38a49578470dd1e8441df1e7a1b3876ab34634b68724830dd2f2d792",
    "smr/app_nuebar_bkg.txt" => "933e99706274f953e593d1430de955e6dc5af74f1d566c3185c900c452de4b57",
    "smr/app_nuebar_sig.txt" => "34976a3cebd25dad3c9d93575a62f7f9e449870c7edb322483104beedf444bb2",
    "smr/app_nue_bkg.txt" => "fce2e9742ad55ba86a24d28bc55d990f5530db89f74fb3fe148b9a601a67c38d",
    "smr/app_nue_sig.txt" => "544778ecdb0e7e851376088be7438eb9a0dad8143a20e081d5283416aa1f2d3e",
    "smr/app_numubar_bkg.txt" => "3afdc37533614907fd89159cf494618c222664e3b4c78f56397ccc6cd64d5791",
    "smr/app_numu_bkg.txt" => "a47d7cc6762bae216f7aff46f0d808721d447fa9bad83e97ec56438a98d41782",
    "smr/app_nutaubar_bkg.txt" => "e4436c5f68521fb6c0856ca841b0a12ee540a2a4c2fc6282378b0a3be1477549",
    "smr/app_nutau_bkg.txt" => "9686f4c54c42bd5e8b22b654d618691d7e21b113a2573fdb028baf034b9be948",
    "smr/dis_aNC_bkg.txt" => "ae6a384bf8fc161aed29726b19a35414e07c53107cd4e8b6aa0f4e9b1a625432",
    "smr/dis_NC_bkg.txt" => "7fdae2f4cd61d6fd3baf138dc9617eadaa5427809aa7b4b5ae0a2a00968ec49e",
    "smr/dis_numubar_sig.txt" => "91d2c3d022e463fdb2f37990a0c73ebfa48b73f9a0abbb538ac3b4fad03a9b35",
    "smr/dis_numu_sig.txt" => "354b74ebfe9c2d60c98d9991e316ca41aa9244eb16ab77302ca7088bbc93b460",
    "smr/dis_nutaubar_bkg.txt" => "652e8a12db70e168a7d428d6efebda8856e8e6ee7ab223338509b4a07f115474",
    "smr/dis_nutau_bkg.txt" => "45fabb987ae6e7a9d2cb8ffac561f4c7cc725d5a1dedb6e924aefe3d45911eca",
    "xsec/xsec_cc.dat" => "bb950e7fb645bc5cdb052d17e2076534be2cbfd24e83ddc851f6d906ece93efd",
    "xsec/xsec_nc.dat" => "613cb9005dae60e028e429f125b537b11bce3f7aaa389d6dc329f4a9ce7f28d9",
)

dune_cache() =
    joinpath(first(DEPOT_PATH), "scratchspaces", "367f8d2e-8d77-4eed-9bd4-73b32ac618fb", "dune")

sha256_of(file) = bytes2hex(open(SHA.sha256, file))

"""
    dune_configuration() -> path

The AEDL file of the DUNE configuration of the Technical Design Report
(arXiv:2103.04797), ready to load: the AEDL files of `examples/DUNE` are
copied to a cache directory and the data tables (fluxes, cross sections,
efficiencies, smearing matrices) are downloaded from the ancillary files of
the paper on arXiv the first time. Each table is checked against its SHA-256.

```julia
dune = Experiment(dune_configuration())   # the DUNE TDR configuration: 40 kt, 1284.9 km, 6.5 + 6.5 years
```
"""
function dune_configuration()
    dir = mkpath(dune_cache())
    examples = joinpath(@__DIR__, "..", "examples", "DUNE")
    for (root, _, files) in walkdir(examples), f in files
        endswith(f, ".glb") || endswith(f, ".inc") || continue
        target = joinpath(dir, relpath(joinpath(root, f), examples))
        mkpath(dirname(target))
        cp(joinpath(root, f), target; force = true)
    end
    for (file, sha) in DUNE_TABLES
        path = joinpath(dir, file)
        isfile(path) && sha256_of(path) == sha && continue
        mkpath(dirname(path))
        url = DUNE_TABLES_URL * "/" * file
        @info "Downloading $url"
        Downloads.download(url, path)
        sha256_of(path) == sha ||
            (rm(path); throw(GLoBESError("$url does not have the expected SHA-256")))
    end
    return joinpath(dir, "DUNE_GLoBES.glb")
end
