# LDA
Presenting a paper we wrote for a student project in 2013. 
The project and the comments in the code are written in French but will soon be translated.

Latent Dirichlet Allocation is a Generative Probabilistic model. It was first desribe by Jordan, A. Ng and, Blei in 2003.
The range of its applicatio is very broad. Let's name a few :
- Topic2Vec model by allowing the description of topics within documents and corpuses trhough vectors.
- topic discovery algorithm
- clustering algorithm 

The paper name LDA_description_implementation_test.pdf is a report written in 2013 by me and three other studends from ENSAE under the guidelines of researcher from the CREST.
The purpose was to demonstrate our understanding of the model Latent Dirichlet allocation, to implement it from scratch in R an to apply it to a given textual corpus. 
More specifically, you'll find :
- a formal description of the model LDA
- a formal description of Variational Inference the training method that we used 
- an application of Variational Inference on a simple Mixture of Gaussian model 
- an application to a specific corpus
- comments on possible improvement of the model to other tasks.

About the scripts : 
The Variational Inference that we used is a two steps optimisation : 
- Given initial alpha and Beta we compute the best variational parameters phi and gamma.
- We optimise alpha and beta in order to maximise a specific log-likeyhood
and so on.
Fonction_LDA.R is the main script that includes the all optimisation algorithm
It recquires the script Fonction_Intermédiaire.R and LDA_loglikelyhood.R

About the input : 
From a textual corpus, we first build a dictionnary that uniquely match a word with an integer.
The input of these scripts is a preprocessed matrix-corpus. Each column corresponds to a document. Each element corresponds to an integer that corresponds to a word according to the dictionnary.
