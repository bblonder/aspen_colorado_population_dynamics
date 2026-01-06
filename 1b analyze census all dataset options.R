PREFIX_TYPE_LIST = c('v1allAspen','v2closest11','v3original11')

for (i in 1:length(PREFIX_TYPE_LIST))
{
  print(i)
  PREFIX_TYPE = PREFIX_TYPE_LIST[i]
  source('0 analyze 2018-2013.R')
  rm(list=setdiff(ls(),"PREFIX_TYPE_LIST"))
}