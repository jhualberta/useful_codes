dataloader = TMVA.DataLoader("dataset")
factory = ROOT.TMVA.Factory("TMVAClassification", fout,
                            ":".join([
                                "!V",
                                "!Silent",
                                "Color",
                                "DrawProgressBar",
                                "Transformations=D;D;D;D;D;D;D;D;D;D",#I;D;P;G,D",
                                "AnalysisType=Classification"]
                                     ))

