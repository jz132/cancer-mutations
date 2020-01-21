table_regulatory_mutation_by_gene <- table_promoter_mutation_by_gene %>%
  full_join(table_enhancer_mutation_by_gene, by = "gene") %>%
  mutate(count = replace_na(count.x, 0) + replace_na(count.y, 0), 
         regulatory_length = replace_na(promoter_length, 0) + replace_na(enhancer_length, 0), 
         expected_count = replace_na(expected_count.x, 0) + replace_na(expected_count.y, 0),
         var_count = replace_na(var_count.x, 0) + replace_na(var_count.y, 0)) %>%
  arrange(desc(count)) %>%
  mutate(p_value = ppois(count, expected_count, lower.tail = F)) %>%
  select(gene, count, regulatory_length, expected_count, var_count, p_value)

table_regulatory_mutation_by_tss <- table_promoter_mutation_by_tss %>%
  full_join(table_enhancer_mutation_by_tss, by = "tss") %>%
  mutate(count = replace_na(count.x, 0) + replace_na(count.y, 0), 
         regulatory_length = replace_na(promoter_length, 0) + replace_na(enhancer_length, 0), 
         expected_count = replace_na(expected_count.x, 0) + replace_na(expected_count.y, 0),
         var_count = replace_na(var_count.x, 0) + replace_na(var_count.y, 0)) %>%
  arrange(desc(count)) %>%
  mutate(p_value = ppois(count, expected_count, lower.tail = F)) %>%
  select(tss, count, regulatory_length, expected_count, var_count, p_value)
