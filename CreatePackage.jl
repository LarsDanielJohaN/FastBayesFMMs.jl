#Author: Lars Daniel Johasson Niño
#Created date: 04/01/2026
#Purpose: Create FLFOSR package. 
#This package follows the instructions advailable at: https://julialang.org/contribute/developing_package/


using PkgTemplates

t = Template(;
           user="LarsDanielJohaN",
           authors=["Lars Daniel Johansson Nino"],
           dir = "C:\\Users\\End User\\Desktop\\ITAM\\Tesis\\FastBayesFMMs\\", 
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )

t("FastBayesFMMs")