# test pipeline to aggregate "same calls" by multiple callers

# double-checking our results using AnnotSV files
subset(ExRes, select = c("SV_start", "SV_end", "SV_type", "ID", "Caller")) -> A

aggregate(A[5], A[-5], unique) -> B

B$ID <- NULL ## Somehow had to include the ID then remove it to make the formatting work...

aggregate(B[4], B[-4], FUN = function(X) paste(unique(X), collapse=", ")) -> common.variants

# checking callers power set distribution
table(as.factor(common.variants$Caller))

rm(A); rm(B)

# plot
sample_size = common.variants %>% group_by(Caller) %>% summarize(num=n())

common.variants %>%
  left_join(sample_size) %>%
  mutate(caller = paste0(Caller, "\n (n=", num, ")")) %>%
  
  ggplot(aes(x = caller, fill=SV_type)) +
  geom_bar() +
  theme_bw() +
  scale_fill_discrete(name = "SV Type") +
  labs(x="SV Caller", y="Count")
  
rm(sample_size)