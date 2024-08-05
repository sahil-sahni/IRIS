
# Updated from original submission 05/22/24 SS

## additional statistics

mat = matrix(data=c(139,160,8,288), nrow=2)
mat = fisher.test(mat, alternative = 't')
mat$p.value
mat$estimate

mat = matrix(data=c(173,126,20,276), nrow=2)
mat = fisher.test(mat, alternative = 't')
mat$p.value
mat$estimate