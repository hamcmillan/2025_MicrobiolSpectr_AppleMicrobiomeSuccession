library(tidyr)
library(dplyr)
library(readxl)
library(openxlsx)

#### ITS gala pulp ####
##### import wald result, tax table, and indicator species results #####
wald <- read.csv(file = "DESeq2_Results_ITS/pulp/all_Wald_results_pulp.csv")
tax <- read.csv(file = "DESeq2_Results_ITS/pulp/tax_table_pulp.csv")
clust <- read.csv(file = "Figures_ITS/clusters_pulp_relabund.csv")
inds_fb <- c("988de518d1836bcb26578673c5a0849c",
             "74c6a5ea0600e02b6e636e617585ed3e",
             "3baaa1a55f67a004e5bdf88baecff544")
inds_pf <- c()
inds_fl1 <- c("ecfed4ed448c412970e5dfeef71342ff",
              "daba98181bcd36fb1590e6f733576f47")
inds_fl2 <- c("02b351f4ae26aba502bf995494d84845")
inds_fd2 <- c("bf00fa4e0f79ee58ba7229a892df5ca5",
              "78baeb4acaf2da94c88c9c8a227cbbdb")
inds_fm1 <- c("8c76b9e28bcc0f156deba52753212886",
              "0c03c7cf289f1f043fbafc14461d9772",
              "1b3eb48f4e7d604c97dbc8f7b9416dde",
              "6b5646add0e0aacad0d377b517427e0f",
              "032c8c83488e8dc1707c235017237ef7",
              "67a311061ea0282e400c59c52375f485",
              "6ce7604e712f963bbd0ae0f6b7a3cf5e")
inds_fm2 <- c("fcf648d3292808e931245b3ef682bb47",
              "c00ccc168f3ec7b9acae23b96aeed71f",
              "a2878a4e252dfa8b1162300ca2040ca5")

##### add taxonomy to DESeq2 results #####
data <- full_join(tax, wald, by = "X")

##### add indicator columns to data table #####
data$Indicator_fb <- ifelse(data$X %in% inds_fb, "Y", "")
data$Indicator_pf <- ifelse(data$X %in% inds_pf, "Y", "")
data$Indicator_fl1 <- ifelse(data$X %in% inds_fl1, "Y", "")
data$Indicator_fl2 <- ifelse(data$X %in% inds_fl2, "Y", "")
data$Indicator_fd2 <- ifelse(data$X %in% inds_fd2, "Y", "")
data$Indicator_fm1 <- ifelse(data$X %in% inds_fm1, "Y", "")
data$Indicator_fm2 <- ifelse(data$X %in% inds_fm2, "Y", "")

##### add k means cluster information #####
clust <- clust[c(2,10)]
colnames(data)[colnames(data) == "X"] <- "otu"
data <- full_join(data, clust, by = "otu")
colnames(data)[colnames(data) == "clust"] <- "KmeansCluster"

##### export new table to file #####
write.csv(data, file = "CompiledData_forSupplement/ITS_gala_pulp_compiled_data.csv")


#### ITS gala skin ####
##### import wald result, tax table, and indicator species results #####
wald <- read.csv(file = "DESeq2_Results_ITS/skin/all_Wald_results_skin.csv")
tax <- read.csv(file = "DESeq2_Results_ITS/skin/tax_table_skin.csv")
clust <- read.csv(file = "Figures_ITS/clusters_skin_relabund.csv")
inds_fb <- c("bae02ed32560e67611d89186aa18803d",
             "4ecefebfee91d1c92b5cbab245111e1d")
inds_pf <- c("2570ec1eb165e1f5a6264548ca2e911b",
             "a0ee0b621db51fe9894fc656adb7d5a8",
             "f3203b95994ca6e1d8121f7ea347dea3",
             "70cebfe7a3d40197bb33a04eadb535b9",
             "97d40dd7111a32700fccec222434ffc3",
             "b6a20dccaff359a64577d6288058120f",
             "5d6b0516dfd5a47e4d6e2be4f5957985",
             "6fcbda34b5f9b3758490edacf6b196f6",
             "327d166553aa3c2273050ca72918c132",
             "05bfde831df4b0303697e3d8c69df357",
             "29a848b6269f511ed03f97dc6da2c60f",
             "94d82a21883867ad71abf46064a63420")
inds_fl1 <- c()
inds_fl2 <- c()
inds_fd2 <- c("0a369a71638f4c13c56f24c4379f16d9")
inds_fm1 <- c("3baaa1a55f67a004e5bdf88baecff544",
              "c00ccc168f3ec7b9acae23b96aeed71f",
              "6a8fb45e374046635b9a1e58508fcb0b",
              "74c6a5ea0600e02b6e636e617585ed3e")
inds_fm2 <- c("20cb2381c7329d2363cd71598db48a88",
              "68888dedcbeabfab48bf3dc093b8be07")

##### add taxonomy to DESeq2 results #####
data <- full_join(tax, wald, by = "X")

##### add indicator columns to data table #####
data$Indicator_fb <- ifelse(data$X %in% inds_fb, "Y", "")
data$Indicator_pf <- ifelse(data$X %in% inds_pf, "Y", "")
data$Indicator_fl1 <- ifelse(data$X %in% inds_fl1, "Y", "")
data$Indicator_fl2 <- ifelse(data$X %in% inds_fl2, "Y", "")
data$Indicator_fd2 <- ifelse(data$X %in% inds_fd2, "Y", "")
data$Indicator_fm1 <- ifelse(data$X %in% inds_fm1, "Y", "")
data$Indicator_fm2 <- ifelse(data$X %in% inds_fm2, "Y", "")

##### add k means cluster information #####
clust <- clust[c(2,10)]
colnames(data)[colnames(data) == "X"] <- "otu"
data <- full_join(data, clust, by = "otu")
colnames(data)[colnames(data) == "clust"] <- "KmeansCluster"

##### export new table to file #####
write.csv(data, file = "CompiledData_forSupplement/ITS_gala_skin_compiled_data.csv")


#### 16S Gala pulp ####
##### import wald result, tax table, and indicator species results #####
wald <- read.csv(file = "DESeq2_Results/pulp/all_Wald_results_pulp.csv")
tax <- read.csv(file = "DESeq2_Results/pulp/tax_table_pulp.csv")
clust <- read.csv(file = "Figures/clusters_pulp_relabund.csv")
inds_fb <- c("b7bd3144ef8660adfb0865788b06eb4b",
             "1517bc324ddc23870f7d5f0d22cabf4c",
             "d7fe05da337d7b4442c6baa2df4ccae2",
             "e0c1d0207a6022a4b81f6a728dcae5ab",
             "aa2a58fc61062baf10d9ae1d28e998e9",
             "4888f12a399b736d232ed0473b16a824",
             "4892d246989a5199a0e4a68164b8da91",
             "d045a33da536c9ed1a68895c710b3c5f",
             "f1d012b83af1853bf87cb11c6279c697",
             "b3cc62a102846d5a3dcaed0571e536d7",
             "fb091f9a2648af2a440d1882fa23e70f",
             "79296f6de2ef4a6a8eef02ef1ba5ba67",
             "48539246f26498c99c5ee72674ef9bb1",
             "1773f8ca97cfdb5eb7680d939ac0c4e5",
             "0306f6e388518a832f37981d61fa13bb",
             "0c2774135fe0c70bf184309639fdbf42",
             "8b8edca944ae3759812e15c726811046",
             "b062c6df0bc5ae1edf4fc91790859991",
             "cdd6484bf04303597ccbe02ca65a6b5b",
             "53962d362402c57d18d7fe53fb9fbc4a",
             "dbc576ef3e30c8d3cd3c6a28fd4d24bf",
             "19a5476f0d54a0eb7cc40638712d2858",
             "6f98bf11afcb6923fff750a581c4d2e1",
             "fc19b56ae03d99436f4e0424019d0eb8",
             "e18a58fc4fcac4385b2052f56ba9c630",
             "e3b8736eb612389a1a7bb2debdfb76c1",
             "bc609506e41f67e3c84c984513005674",
             "18128d2e28bc51db581d003f8312cbd2",
             "baaf68af5d9bd5210cbd466536ee3f30",
             "4b94439bb6f99d414e1230ddc444a935",
             "3f28a7a0ef8a28c1c753a3c1a134d31c")
inds_pf <- c("c429b65ad7b95849da841c9a4239f06e",
             "208368417cf91026ee691630fef5252c",
             "bc876e7885cd5839dfc23670f898fb2a",
             "502e8a712e3ae5494ded27e6482cea55",
             "5da0c348da46f51f3ec93c4d590f2223",
             "07feb6e1850dc105154036186fcde3db",
             "bc1b84288e39c73b2d8356d391c9a3fa",
             "b760288117bd49ca72bb8c8f3717bf79",
             "e39ca4dd939d60e45603471795f698fc",
             "b9b0f4ebf187727476e52033bc6fcbc5",
             "e1c32e73b241c1aa92010422674bf406",
             "ae63ba67aae0c2c42e5f051f9fb1d7ba",
             "0898d667077cd396c2514d10b0c5c023",
             "20c2397fca0fd722c0446f9805837380",
             "2bc3fbb77968c4cdaa13b868e372b439",
             "ad320e465c374a5c083772457cafb0f0",
             "c934e224ba4422f2f37388d30e7c195a",
             "bef983dc450ab39fdb44b10c691a7902",
             "9155b1b596b2db9ba8d62259a424d65e",
             "5e65b0e244c9ea4bd1137842b28fba27",
             "96dcf6a010c2927628255277241ccc2c",
             "b196e34e5611285a5243d0ff99ff26f6",
             "b3d76aaef27fdead0bd60b019d89662e",
             "9292d039fa55279cf6336b7f7cdd542d",
             "4411bd3d081f190478d60a8fb06ce7e8",
             "a1622c7f144e462bd3edc0ea5fc64ca3",
             "787b32f4b5d15a6037431afab3a00036",
             "9d6c33ec04f0c27ed93bcd8defd0e051",
             "1e2178df902293fe72781ab74e5551be",
             "c3d6c4ee62c1543fe5769f0f06f28f68",
             "821e1a9200f8e3fd31901805cf673a19",
             "62869493b0eae6b92153a75c1e9d519b",
             "8fe9ed11e435c700f9152edd4213aacb",
             "15b19285f9a755aef07076c6f18c701a",
             "261282d3f48e3490b9fccedd1dc7dd9c",
             "354559b2fd0b1442a9a07b48fcd620ba",
             "7c616643f2bc933ac56205af7c0f60b8",
             "e3f97c22c903201f99cf44fcfe790370",
             "919d7bd67e8f36ec07c3ada12a4e0352",
             "e3128fc690ecce7bdc8c503102e9f244",
             "68ca9743ab6ed940b58e9bee2d615eaf",
             "c540a631da77d7e2f5efe7a8aada2f16",
             "841d45d1db558a8909f142100597a7a0",
             "35ae5f93ab0d017302b479ca8fb0e618",
             "5d425f42e9d2ecfb212f860b2120315c",
             "bafbb060659c1ac720a5212d5c1638f8",
             "33782dbf567b1700759a181a9464d3f8",
             "9107e526ae87f403523dac3c5580f180",
             "c5f90bacdb69d585bc439ba2185221aa",
             "0644404075059c60754f3c4cb9132069",
             "c061600b3c9fbd2e015c51c4cf9c217e",
             "aba86589778e4df4d6b9355a51c6ea7d",
             "696150333b2af9ca53f97fadc6c3b582",
             "888e002f3f68b91b28040c04612e976d",
             "5091c329ae0f70c8d7e2dfdb06c4d2cb",
             "18785614a70405ac138ecc8fe05e7af6",
             "10e0329dc70656d2034d72eca9c1c197",
             "ee2373fcb553a5f083801e17825653ec",
             "1000084d9c55a258e7a84121c55abdf4",
             "9294a9208de7ae50d655f23fdaf1fe79",
             "74d7a898f612dd08825367c24e411fa2",
             "b0a06ab806294afbc4b1d0464d235eb3",
             "5f2a36bf7fe14eff16f75a98691fdcdc",
             "ed2d129013d67899562368623f37d414",
             "31ee31503b8153586f9c9165caf5d458",
             "f5f02dd19a37eea1275dede4900bcb6c",
             "f8656365210009b995d2994ea9e21627",
             "93abb9cc764234a274808b080050d361",
             "fcc473c9cdbab45aef52525a28135374",
             "34612c5fd4d9090d3bfacdb11a9e3f7b",
             "3f42634c8cdb041068725ce96a4dcbc1",
             "53522396782b1315911888df672d91aa",
             "051787a54755b4b148df20ee5da36b6f",
             "dc09dfd2cbada12c3acaee3a656da795",
             "6fb57b3f591c06ef265288927280204d",
             "8d531059671fbc0987ff2550ab78ce7b",
             "38e93c760110ae1554c2c758aac30750",
             "42d018b0c8aca42bdf22def2d3132e9f",
             "82bc923d4c805609eacaf2bd0b683d94",
             "b257f68c119d2a02e709fca0acbac0d0",
             "285522c774bf6d97e6424f95a048f6cd",
             "bbeec0fef91f3f3bb95a451ba7ac4bc2",
             "40159c0d287968608e0f3bfb52aabd7d",
             "ab2228c285e58143e3937619a17ca59f",
             "c2bf61d55dcbd06677fdf9a6db500a74",
             "ff0ff1fc050172a114d1a4abe49e48d6",
             "fa0840ce3b7a30779b7d063ef228b584",
             "a4697b7b09186f1afcc676e2037d734f",
             "f07818275d231b8bcbc85d12a62a446e",
             "6da60c5202f851769494c809031de291",
             "9994921cf9bbb84f0dd12e7fe67f7ea8",
             "993e798ce046afb7988b81752d7a2601",
             "fb87245876c4601269e6bafd34e25476",
             "d02fa8e828d8ecb93580747f771d5efe",
             "e5e6244cb2bf5467137663441aa52a23",
             "6f4864fe4530522cb5bbe3549b7ffbc8",
             "4ae5e913d72d8ac513ef4b84ce2062de",
             "75e102dd2cef5ea7c8ac4f6e954b3e19",
             "89265b066bdd5c00f873ded76d503ed9",
             "a2a54507fdf0032b592548f4a537c4af",
             "f55053d91be4dd9bfc1e084afbeeb3d2",
             "8ba325502172a82fe69139c092a41f30",
             "a812ec4aaaf745be142134bc9e6ff6f8",
             "3dc9cb3b4f00ce911d4b447446873b39",
             "8041305f0579fbcb96b79686a076a676",
             "1adf3a8a77fa70c7746bd69d9118120c",
             "d71b1b41858c80da61fa9f049d38b68c",
             "57a58b4e2a8f36e7080138b28a7840ea",
             "8ed3f7672a93c0971e1a481d2f806a6a",
             "aaa08012744c21876e1c0f9e9da479e9",
             "f8b6332c392e62acbd674dac63337a36",
             "4ebd718b6b5b616aa8778741a47c0dbf",
             "4a889494166eed75e234b179d37a1047",
             "3d50ffab3d2b65dfdfc89264bb95c32d",
             "f3de0a1dad7fd10675b116ef59bc4f2b",
             "0c2f6f4c04c0f2ce8e0cd478645fe1ff",
             "727c68dd7746e1db7255c74ac6987018",
             "09c49e03db0d8948b8b01ebe3507f7ca",
             "dcd9163b0b5d26f7f3e99017805ab23a",
             "facaad7b2f0d0dca1a40b0ce10c37a45",
             "a2a65572f398540e2deb9f671f6f6100",
             "abaf72d716d95cd9bcaa76cdbd1be21e",
             "cdcd0b4c0d5188647d43393b7d65b88d",
             "43dc75d6424595b95d39acd8521cbc64",
             "97386bce4f48091eb27a14fcf10b6885",
             "f1de9f421542bf342680a8a1eb50a226",
             "c321c8f0e604096e3e581cda31ac12b8",
             "03c1bf5ef6df57fe9a93f7f7023e5bf1",
             "43bdc9c07e7747c762d6d185d038065f",
             "1849d6fe64e4f72d94e02b33b0e9f125",
             "0ee154cb2b34477e6b91e1c26086c5f7",
             "975625c24ff69ac6bfd51f1a0636fcf3",
             "42a5164343e88f24e70f4d93bd662461",
             "c073ec2dd0eb4b979aca7b7addc5171a",
             "c276461d20ca3fe0141be748d0e745cf",
             "28a2784faae96dd80471a4b174af6fd8",
             "d4ff8438672d8a68d355867ff8a01499",
             "aebe511585a1097dd61e78ac56c729dc",
             "ab046e73baa2a26c2d70758979ce0485",
             "b640f2f8756b1d26fbaaa1f3053e2daa",
             "f748a808729ca80d64a4e6cc48ab48e4",
             "06e0f05758c9f55703487becc6ebc352",
             "5a97f4c2bee806cca2802fbb141e5d5e",
             "d3a06d59281c5102213df4d003a23264",
             "359462dddaffbb33b5e91144c77ddf6e",
             "65ce1dd734930ff213fe2cb8bfc14f94",
             "8939ca7e158482cc7005001776531c22",
             "bfca4d92808583e17e6f2ef286595229",
             "f97e4c7b9d9e673ce32f8a02523b628d",
             "2b7f7a2cfdcfe0a68b5ed84e246d087d")
inds_fl1 <- c("758f71e9312cbd3d34b5cd77801eac2e",
              "ddca2aad3b89700e48f170ac2bbade8b",
              "1281a215c0c390131cac85ad1a9cddc4",
              "b8b74a6606d4202a5dc86d6eca12d3dc")
inds_fl2 <- c("3b226506634113f101aeced89ab97917",
              "7610e232e127c76e8358460411f3aae6",
              "186e8e86f025a6bb5375c1f9641efe41")
inds_fd2 <- c("32b44e941fe682173c01f3d1cc72146b",
              "d07dd984f103e7586c7fe7c027a50067",
              "4f4b8775a66491d31c87f5b704162142",
              "a87b99db76183dd7ab855a9ec34e26a6",
              "d64baa7e97cb499e1ac80dcb2fdbf1ac",
              "e39dbaefce3774cf3c6f38849fb2bcfd",
              "1464cb49af73be47ec803262bd9623fc",
              "2953367c5a22718206854d74101e4e35",
              "b478431dabab00865efcbcc93161b463",
              "0a92f2eef1a2a18504fec06fa76ef059",
              "05fe3b6c8a5579b1b0adc54b867bb4b7",
              "b5490721c6336d7863a24bc78d9b6d9f",
              "888951e7bde7438da05d26120b769c08",
              "40a4f0fe0d5a48714e975721d0573926")
inds_fm1 <- c("8497e2776fc351202a443e2a93dd8610",
              "45d87b9416b1a5573bd6f1c7d764a57d")
inds_fm2 <- c("9870ff2232ecd0ca58a8504babd0f38c",
              "aaa7282c32d99183368853ba784c7550")

##### add taxonomy to DESeq2 results #####
data <- full_join(tax, wald, by = "X")

##### add indicator columns to data table #####
data$Indicator_fb <- ifelse(data$X %in% inds_fb, "Y", "")
data$Indicator_pf <- ifelse(data$X %in% inds_pf, "Y", "")
data$Indicator_fl1 <- ifelse(data$X %in% inds_fl1, "Y", "")
data$Indicator_fl2 <- ifelse(data$X %in% inds_fl2, "Y", "")
data$Indicator_fd2 <- ifelse(data$X %in% inds_fd2, "Y", "")
data$Indicator_fm1 <- ifelse(data$X %in% inds_fm1, "Y", "")
data$Indicator_fm2 <- ifelse(data$X %in% inds_fm2, "Y", "")

##### add k means cluster information #####
clust <- clust[c(2,10)]
colnames(data)[colnames(data) == "X"] <- "otu"
data <- full_join(data, clust, by = "otu")
colnames(data)[colnames(data) == "clust"] <- "KmeansCluster"

##### export new table to file #####
write.csv(data, file = "CompiledData_forSupplement/16S_gala_pulp_compiled_data.csv")


#### 16S Gala skin ####
##### import wald result, tax table, and indicator species results #####
wald <- read.csv(file = "DESeq2_Results/skin/all_Wald_results_skin.csv")
tax <- read.csv(file = "DESeq2_Results/skin/tax_table_skin.csv")
clust <- read.csv(file = "Figures/clusters_skin_relabund.csv")
inds_fb <- c("3b6dfcacf16a90d9245f7a51b8413a54",
             "e487b5270823c679da1fe23768f0eb3f",
             "faf3fb9e0acc8bb6ce2f2830bdd6e516",
             "dfc3599b82444dedb43ed0686dcceb18",
             "d138ef02e573e5fc0a63b69f323781ca",
             "8dfc332291e678557a7a9b75570d94e0",
             "a786abfe6cbf7835bf5b96027101429f",
             "88ff88a3a4e43b7120c3529e8f4674b0",
             "b023e3f22452961aa8df340abfdd89cd",
             "781a94320f21af428fb308a8e1386c91",
             "69898b90eaa01bacdb31346f9bbcf160",
             "379d83d55abfb88f01decba36fdf2be0",
             "7fc9c04ff88a2a72a739f51b1784bb60",
             "2f2a1e9f86ec2e10bf122d6d55b858fe",
             "92d220402fbf1de449d87bb77564c888",
             "e9d71e19d594780b3483fd497d67904d",
             "6c3172193eaf4e7a67e931b490180368",
             "e992ccc75b16fbdecabd73d3d08da58a",
             "779fe775c26124d8646bf50c697b88fe",
             "f8bf9e655944435a0a0c3857ce81d43c",
             "5cf59daf13dec7a37ef537e16d027c9f",
             "dbcd1c5d1b1f002b3b84e49f9f9ac4ab",
             "7e70d2eaeb752d87f2d14936c15d695a",
             "00f49789b6f6ea817e7e771936bf60be",
             "7a57568808593785b4401f815077fc63",
             "e9fbc5638db734850f83cc3b2d9d7fd0",
             "6ce3a8d67fab655379e01bfb394998cc",
             "b8391290f14b3bb5e179f4f63aba91e4",
             "3a0f9cb485098bc6a7ad1244ad57b313",
             "9f003b52559f2fdfd88cdb2e0bcba2c9",
             "52319a5c30bdf69523d2aa6590690c33",
             "5b6c12df1d137bfa6e5c2bfc5512d7e1",
             "184e17ad2209c98a487e2d7ea34c2d85",
             "100e8de25771e0ceb8b96704648d0dd8",
             "06dc028f61252624c892dcb62eb6e5a7",
             "4e350caf3af4d7516fb2d7ed406759fb",
             "91b521f69ec7f52128b679b3b4c349e8",
             "90b6c9e9aebad202c558e771f04978ca",
             "c061600b3c9fbd2e015c51c4cf9c217e",
             "16f9b022a613318b14b33b836c482d6b",
             "6e967f06fe5a32baf133706d3fc04e03",
             "7e65eb8a15cd82c48c41004f0f1ef17a",
             "58a73cf13172d8d35ec4aba726f90d0a",
             "f9b701365f07415906f51e80a22c220d",
             "a8fc3af81601116e217c0cc5ea9d43c2",
             "aaa4d56b62b50f5c632647e0819327a1",
             "1f54a21baaa23e6076dc8d27c021e402",
             "c83dd3bf19c24dbf825b33f6f14947a3",
             "cdcd0b4c0d5188647d43393b7d65b88d",
             "e3b8736eb612389a1a7bb2debdfb76c1",
             "55279e1b096ee5c3459f1628134fbb14",
             "2d873edd98f77e3ddad442054a921bed",
             "0a7cf7e2c2e487165f70c060e312a820",
             "2c11b889664708d1e377ed01740f40e5",
             "9dca40ec15607f2e3c536c12a9be6499",
             "426e70b3aef8baca5f157f89d0db83b3",
             "bfcbd21fd7d231c1dbb4a244486f76bd",
             "957b021f727f2065c7cdbc08646bebff",
             "3cdfd20e26b4b29e196e8e3c6411943c",
             "68ccbdee013fcf932192f31c18124c65",
             "dc24422870faadd3a12277a02e1b042d",
             "385d9bbfc2844af5a2eb216dced40a70",
             "c692c855775e2c05eb2ed38382ee1bd0",
             "a0fa5bea09eb790622df9bd185137893",
             "917c0649e005d4c4ac23cf915dec24a5",
             "27f09a775b0416a33eddd6e0adaef5ee",
             "6b57f261038b434bad79b0d6c5ebb083",
             "0d454d8d596049d6746188510efe014d",
             "2349a041c661ede36101305b24019821",
             "dd0ac2fb5cadbe2499b83dc3d203bca6",
             "9f874970f1e27001a56fc7df2f65ac9a",
             "8656db6865d61d94102711faf9792bc7",
             "ab046e73baa2a26c2d70758979ce0485",
             "ddefce20035fdb4b2ed12884c865d8c2",
             "4c67c18f59fb5a6c711c58de56f430ac",
             "220a07194985faf86da06bd0bcabf721",
             "ee85428b4a2121e3bec628668a9433a6",
             "4dca5bc8b0b91e806d7dad413bf49eb2",
             "133a4e6e3bb42a4175993821aa7fc0e4",
             "2462bb13bd615701cf5d0b79a0a69969",
             "4ae6a3a9e1007cc08ed573ba8eda0a28",
             "0a8c636835d37440c8129d2b3321bcd6",
             "caaeec287a94051c3a32d2bc6ea9e445",
             "63ee3d5713ead82320bf124582334273",
             "702fd3ed00d132fb101b722cb0a62f04",
             "7ad7c1618666ce678eb9b30d6fa2ef3e",
             "939db6fa73104893b4bed18b30f20d06",
             "a8edbe835fb45fe2757afc2c798aacb6",
             "3612fd75c4163949d62e0c86df909c4b",
             "8d0ccbd006e7f0804ccd2e344b6639e4",
             "1a453cfb9a2376460453056acc34d0d7",
             "0b97876e4f841af80ee9f34df1849963",
             "e58e38dafd0561e92e1756c7b0e52b67",
             "921417ff5d86eeb0a341914e9067456e",
             "4b94439bb6f99d414e1230ddc444a935",
             "0d1cd25af9afbc9dbe4dfb92ec44cf73",
             "78f2c303c6eac7d4d70541eaa04ec51f",
             "a0d604e8cbc18f96d0d0dc7fd8f273a0",
             "0da48b990716d8a7e37f3d2fb2c83932",
             "c582af7a6b70529d22616aaefc4fdef8",
             "1e9e889155608906c3da26a11401583c",
             "45048842071376155f523f3a2cf1d7b1",
             "1fd99f7353bc1f456eca47fdf71cf472",
             "3a56e39ac7fb5066f491125d4d7cf569",
             "b44d70bad79635acbf1413be6157df6d",
             "1b83623c596132c0594576daba38c87e",
             "e46470df7142140a7b06b7719bdae15a",
             "e0774561c49ca646bee06699de0b26d7",
             "ffffa274204969081fe0852a639f59ec",
             "9e413020167322edf55fcc29d83d6031",
             "275b208f94c68d3c0464daac87d6baea",
             "2318e72ba025af2353a7130517d6f6f0",
             "cf6685808bfabe483dffe143b63f6130",
             "5433cfc3a09ab266082fd51a152e3079",
             "92a57cf140d39ec5f2451b1fe25f0ab1",
             "37b31b4e84e1d04e151ecf982e164399",
             "d36c4666632b911c436bf5067b50cabd",
             "6dcaa851a330524e58f566db9eaf2f2d",
             "cef5d05f19ea98a20ee0aeb03dc6523d",
             "8fe9ed11e435c700f9152edd4213aacb",
             "4cadd12c99d91fad864d4e0a32d6b3b9",
             "e0b7064fe8d7f4c46e5858af838afe02",
             "77e57f76b335d4a5278debea37540b07",
             "5c3c4b73c6be148fce75f2fcfc20954a",
             "fa406f005130767eb95ab816183e032e",
             "33f4899b97c4601009136f9f2cc76608",
             "7610e232e127c76e8358460411f3aae6",
             "d7fe05da337d7b4442c6baa2df4ccae2",
             "7e4b29a6347040516896f34ea568d83f",
             "4a889494166eed75e234b179d37a1047",
             "7f4745422e164874bdbad37951736ebc",
             "a1b573f7553e1efb7aa5ce3202146b10",
             "2cd4af056087a86013651dfc94e3eb22",
             "68340df9b5c274c9592ed80250dbf41e",
             "2b269d2df875e13ca1a039abcf237276",
             "6ffc79677ce44a751c5297f39a4d12be",
             "3641d63f50782a14043090fcc6ea56e2",
             "7d99a792a301a8d4b3fe0d9b07117ce0",
             "a82f65ba06f79f14134f04fde8305ee3",
             "5c37631feea07fede9497146b54b3423",
             "5c706d2a62d588d59a1353fe4b011e24",
             "dd01291bf2fc3681df2f24e57d5192d9",
             "7f63e8f9435f310148a8acea1d4e2d69",
             "33c633501504f711c7aa20f702da2d51",
             "c528e0f77ca6fc6f066c6c0da2114eba",
             "791dfb902ce2ff9183223b8e12db92b1",
             "75cdf1986c7d2efa61ccc96d1f388eb2",
             "44280e82e0eaa6621928f71fb891b689",
             "9686449e97d29282d1e62a219f343e77",
             "f29741097530f786f2b49b147bd6d304",
             "06e0f05758c9f55703487becc6ebc352",
             "0fc901d7fa211ccc041a7b4f607c71f3",
             "d808c5a46fdc72b50cc5cc1e00222ce6",
             "1773f8ca97cfdb5eb7680d939ac0c4e5",
             "a168342137ae3e99037d0d562b2e6b55",
             "71522dc1ebbb1732d00957175ec4ed39",
             "4f3f5bd296bb10053e7363d252555eb0",
             "4c2f69b3ccf1f493ca856552a49c4fe6",
             "35585716ad8df1b7d7fb49d2b7dbd8cd",
             "cdb112d746a4c284423b063c9f2e3a17",
             "0f8e831e3daee995d422554ce9846104",
             "d94dc8c6c74f500c1b31868634c4e3c1",
             "2f3b51c06a3fb48716343cde5fb19681",
             "07ea517698d4a013d462538bf14c6b27",
             "8b83a05063a8d6604518a8ce1462538a",
             "412b37f96455d0c3deb66433f4fb7e50",
             "526702646fcd7ca79ba495e4fa9ab3b6",
             "e42d4e0921d5ba1090fd409ffc9c2abf",
             "61f798777f8efa602862411f9f318985",
             "32c7450f224ab6cbe4cf8096cbe885a1",
             "e078a35b9c09a6247396e1d43e5dbe0e",
             "8b56a6920253643f65e16554bf57f67a",
             "f45d4fbe4674915459d2c8b6e12a7607",
             "d8fb9518f80f9ca4160ba9b89bc3865e",
             "5e4421a9fae579fbcb709f3db63bf5b0",
             "6bc0ff70dec058dfe3aaf121059fdde4",
             "f3652dbfeb53c929f76696d43b2cb04c",
             "6e066d787571cb3268083780220eee7d",
             "9705bcbf5c803b063a9bd44d2507c6d1",
             "adfdbb568f835c1ce6f26aa4e8d7784c",
             "0ae5d25af84012a7ee785ef7d1c660de",
             "f77ffba3cae0681059458aaea5bbd1f4",
             "fb32dc6924fe321401aa2a0dfe86baab",
             "078900a7eca13f15d09f487615e824d9",
             "99782763fdc0df9e1b8309b14ed9f8f7",
             "920f1bfbbe6c3148616a99990d5ce9c2",
             "0647d312fc1029033522448422f5f7f2",
             "8f14982668982c6dc1bca9f8aa24a043",
             "b89614baaabf5596b49c2ca6372967f5",
             "989e38c1c565434e8cfbb4218ca423ef",
             "3b3598990d43648c618674eb8eef30ca",
             "b1e09fea141341ff4ec31e3d3373e717",
             "bfd92ee9b828937528779022173ed619",
             "0fb4b2158ceaefa574b4b8b58a4709cb",
             "dddacd36157174ad5ca06eb74d0cf12b",
             "70912149d418e5ba09bcea938697263b",
             "6fece7bb9ea748419ef041d75f0f7169",
             "9710836d70a16c2a1d486f43582c8b82",
             "31829cacb16a81c69a48e7f134d08d5a",
             "44d5c2d82e58dbfde65a12139f8c569b",
             "aa9af71f8dec405f259fd2cd9b494576",
             "222391ae7968223faadf714224998595",
             "d19031c88bbdc2ddd60e5b2ddabf858a",
             "2fa3269ab45163f7bcf4d4b33f05653e",
             "adbbe833ef26697bf2f56f6383a916c9",
             "9a255bec7bf2ec99ff11c14ab94a09c4",
             "15b4de9816bef0b4f582ee708b36e9ee",
             "a18fab2c61fc24ca99444f775223a9bb",
             "5aea4d5196ad7f981169a347df44720a",
             "d8360e1e9c56ad84b89478853ee0a29e",
             "e29ce760fd94ddc35fe372330c5129d2",
             "224c2ac4f07f6502fda847f307b77462",
             "e557aee01b6d7c522507c6d81572dcae",
             "13405ccb69da824e2a156045b5cf1889",
             "e7f383ea45086c1107498a7dda818d55",
             "5cd278c1e93981aa83ba9e9d3fdd01ba",
             "3cbd7192db9551b67b9d347a1c5fd7b0",
             "c6507d9ae0cf04107d122914f6d2c3d4",
             "01a7659811229537bb9f6d5983999893",
             "66319113f7f44c99f589a848fc699a6c",
             "37006ddd24e1830bdd6d385b26282cc4",
             "513ee99dea092d8a95a537c09290cb7e",
             "58dba037c843b6bf30cbf6a0c8edb534",
             "b18f33d969ed3c8a2b02a65a6dff33b1",
             "c87eaac2daabbb5b0b8307ad22946f01",
             "9a840c39f18e4fb5070dabe281e9575b",
             "39de912b3b8f17a0bba3b56a9ac44f07",
             "464d595bdc30aa72d2add3520978e1a8",
             "97c83b4addd51e1ee72305caeeac1554",
             "4fd1c187978faf5dfb483eb268f6d0de",
             "c4507f41e59c82939f449d644d21e316",
             "40159c0d287968608e0f3bfb52aabd7d",
             "4cda443b18a955890a2dad5f27fc3186",
             "7c6809803718c8e56b973a25e6525b76",
             "ca7dbe2902301b22dd3cda3eeeb31123",
             "362a6f3e32c29b24bd22832417a0424c",
             "da12fe5216d7f07f56e0c54629e638c4",
             "4cf03ac5041b8d4184246d31cf6eb480",
             "42b047612f05822f1f78d0d449af0698",
             "5f1bf26f281040c2fe1b763be3c5cccb",
             "6f68e13aab22072af56b4dc4cff4d51e",
             "1493f216b84a4af7a66a07a3c948d67d",
             "5a1403bb1b558e73f62d24dd95b35c5e",
             "e54a3fb5a7816905598101e0e490a3fd",
             "7953adf7d57233449f269b08955cee71",
             "05df31a510c915bf63fc7407162e78e5",
             "70ca30fdf8e83f5e7fbc48f207ea55a0",
             "3180719bb859948ceea7c625bbd4755e",
             "3dfdfee079d0d9f3e9f486261f17a156",
             "b2477b592c18c5e3fde4d20e8c0c2336",
             "767681eee70c6555ec2833905af3046b",
             "535bcc31b965a0e416317992c3cdcc98",
             "d1b7fceca67199708b448e26d79d890d",
             "a25a401b82198f73faf5c28c2be5542a",
             "d2596be6866230505308f2b9e47574a3",
             "242df93101f8c77d678133ba0807ff4f",
             "78344b974459c58bfe4a037d05bf40e9",
             "ed7abcd3a268861c4b05497a5046e96a",
             "e251e612658ee9c209b972bfe19ea79e",
             "4a8473cbaa50db75e8a8caaaac3a36b6",
             "8efd8f7deb873dde4da484f46ad839a4",
             "f2057f104d5989c009b7a40c5c323f4f",
             "dc52055f555e47cfa5e4d04929b3dd8c",
             "d85dd3f5e3b48bf8e8d50f499db60e62",
             "2b7f7a2cfdcfe0a68b5ed84e246d087d",
             "75ad55fe6dca7992c8832ddb3fe8cd2c",
             "5e1a3ef51ae7bbc168c5d00ed21760d7",
             "5a8b6093ab9539a0f3697b496f01bcd0",
             "055c99ddc6fb825b1ecf6accd44360a0",
             "dc3230c35b6a52f44bb3c1fe6abdb864",
             "6d166847c9425d3744828e77b15d2443",
             "ab6baa73ed571e84e352139134c2a5ec",
             "600e3af830326b15f4cc5db35a1fe8b2",
             "6fed37f0fa301e02c61329af870681f3",
             "80d219b87241d6771b10db8540038ebf",
             "525bb1c1c79aa4cbbf4c2869fff61130",
             "10fb3174720b79f71f138f136049f94a",
             "7fe413e5d26558b98517c82c24edb1c0",
             "0c4616b686013ff576725e337f90d09f",
             "5463fb1967055a54f80c174ddc903cdf",
             "696150333b2af9ca53f97fadc6c3b582",
             "8680280f7a3a762fe6542986ce88561b",
             "95b254dbd95168a1422916d7eb353a29",
             "2d225fd28e6d013eed3906003234b0b4",
             "a2fa0290c9d8f19d33c3bdc6d07b3755",
             "98a3644131982eef5ed37de995bae279",
             "ff89ce90c1f4f1322abb08941c9bc68c",
             "db0ed4c3659d19a822908722b885a2b9",
             "179429ac7b6506d0ae1e6c3d2d5ed1aa",
             "6082750f89f897051cdd522a5bffdf40",
             "e5dcecbedba36b8884adcd13ff430a36",
             "41040c9994b2aadaf2b2f1ec6956ed67",
             "a69cb8ded070574a2a995eedc4af96b6",
             "9959a69cd704084a9cb7b361690b830d",
             "7e65d22de3a30f67b79acd518a89f85e",
             "0306f6e388518a832f37981d61fa13bb",
             "1ff510ba917e3e58e09df7f7e41c931c",
             "7db9f340a929d69937bebfcc0db7e38e",
             "3dd9d9b687ad7da6fd77f63ab4d85d1e",
             "d55351f72dc893502aa046f2a78d5064",
             "5d34cef462b0b582b8c85e9df016124e",
             "4a9d33f07a86ea4477048d84148cab41",
             "43bdc9c07e7747c762d6d185d038065f",
             "b20168c00fa145dfeb96b3f0308a4247",
             "b58cea0e23560ad8fc434d614322b045",
             "f9c7f9bf285b6c4ae1b396dc22f37651",
             "911b0f4ff8ee18d37c5b5ab7f6b6fecb",
             "2ba779cd9f6d7c6e35196223ee27c676",
             "648d6afd42328d617918b9d107752c16",
             "22847d7156306b266dd1c65c89064184",
             "604d64d23758cfcddfc1d3e51521f4c1",
             "ad0aa15c4e39b5fe0d84d6c316b7dbbb",
             "2b39ea58f8ca2b7130b5a8afd6af4694",
             "7d42406782e9e6ad54f3cde48c191e4e",
             "c8a4f5490bca7869656fd7259f66f509",
             "ae23cc0aa6c7f0f4862dcc7f361b4b87",
             "0dee8d49753e8d5b1bdfb24efe017ebc",
             "40e2ac4a8c01d122f6a6b85b13cd3114",
             "c5ce1532413bd9eeb7476bdd8fa04f3a",
             "25281358b0a6a69d45d78f83c29a5e42",
             "7e80294ea011e764dc35ce58fd4bc6c3",
             "1517bc324ddc23870f7d5f0d22cabf4c",
             "dca904d7e9f72f2955a6d501dacadf6f",
             "28cd0c25afff520448ac95c6772b3c74",
             "228ef752c263b6892cba1497a6a37f9a",
             "26ddfa2084097997fdc670016cb18adc",
             "d798ded5e9986c9b6fa85adffa524837",
             "ef6e15bde39ec21beaceef72018a82c5",
             "3be6c0e27f1b26766af115544856559d",
             "951b0ade33fa2a5525cbb0d9a340045f",
             "d9b89cd0edaca2da25c1ca19e5cb0275",
             "cf6464a9cd44b5d535c815e3ec222471",
             "6236f2e3233a240baecf87e0033bf03e")
inds_pf <- c("62869493b0eae6b92153a75c1e9d519b")
inds_fl1 <- c("ddca2aad3b89700e48f170ac2bbade8b",
              "83b99fcbe972a1df2db695f6733072b6")
inds_fl2 <- c()
inds_fd2 <- c("dd1a516c11812be6198001ca9a80b786",
              "aaa7282c32d99183368853ba784c7550",
              "4fd5ba1baa0f617f1c3f44a783c37c40",
              "c4e449eccd82d5990b6fc127cab7f758")
inds_fm1 <- c("e6ee8857a3d1e436ccb0def13492f699",
              "2a0ddae6cec30b6ad8e742143c8b68b8",
              "e05a6dab445d2978152490ddb10666e9")
inds_fm2 <- c("ea63659fcfc6ef973480f458590cd0cd",
              "fac1450217bddb23f41304aa255ac609",
              "67e003eabae5d20f67d30bc7b9692113",
              "29593b8679ce9413b0fbb2e1b272ea8a",
              "143bd0e551e13b953b27da731c8ab7bb",
              "1540631d68a1cd972a5a7443ceb373ca")

##### add taxonomy to DESeq2 results #####
data <- full_join(tax, wald, by = "X")

##### add indicator columns to data table #####
data$Indicator_fb <- ifelse(data$X %in% inds_fb, "Y", "")
data$Indicator_pf <- ifelse(data$X %in% inds_pf, "Y", "")
data$Indicator_fl1 <- ifelse(data$X %in% inds_fl1, "Y", "")
data$Indicator_fl2 <- ifelse(data$X %in% inds_fl2, "Y", "")
data$Indicator_fd2 <- ifelse(data$X %in% inds_fd2, "Y", "")
data$Indicator_fm1 <- ifelse(data$X %in% inds_fm1, "Y", "")
data$Indicator_fm2 <- ifelse(data$X %in% inds_fm2, "Y", "")

##### add k means cluster information #####
clust <- clust[c(2,10)]
colnames(data)[colnames(data) == "X"] <- "otu"
data <- full_join(data, clust, by = "otu")
colnames(data)[colnames(data) == "clust"] <- "KmeansCluster"

##### export new table to file #####
write.csv(data, file = "CompiledData_forSupplement/16S_gala_skin_compiled_data.csv")


#### combine files into one excel ####
gala_pulp_ITS <- read.csv(file = "CompiledData_forSupplement/ITS_gala_pulp_compiled_data.csv")
gala_skin_ITS <- read.csv(file = "CompiledData_forSupplement/ITS_gala_skin_compiled_data.csv")
gala_pulp_16S <- read.csv(file = "CompiledData_forSupplement/16S_gala_pulp_compiled_data.csv")
gala_skin_16S <- read.csv(file = "CompiledData_forSupplement/16S_gala_skin_compiled_data.csv")

comb <- createWorkbook()
# Add first data frame as a sheet
addWorksheet(comb, "gala_pulp_ITS")
writeData(comb, sheet = "gala_pulp_ITS", x = gala_pulp_ITS)

# Add second data frame as a sheet
addWorksheet(comb, "gala_skin_ITS")
writeData(comb, sheet = "gala_skin_ITS", x = gala_skin_ITS)

# Add third data frame as a sheet
addWorksheet(comb, "gala_pulp_16S")
writeData(comb, sheet = "gala_pulp_16S", x = gala_pulp_16S)

# Add fourth data frame as a sheet
addWorksheet(comb, "gala_skin_16S")
writeData(comb, sheet = "gala_skin_16S", x = gala_skin_16S)

# Save the workbook to a file
saveWorkbook(comb, "SupplementalTable.xlsx", overwrite = TRUE)