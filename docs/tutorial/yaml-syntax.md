# YAML Syntax Primer

The **3t-seq pipeline** uses the **YAML** (YAML Ain't Markup Language) format for its configuration files. YAML is designed to be human-readable, but it has strict rules regarding indentation and structure.

---

## 1. The Core Rule: Indentation

Unlike many other formats, **YAML uses indentation to define structure**.

- **Always use spaces** (usually 2 or 4).
- **Never use tabs**.
- Elements at the same indentation level belong to the same parent block.

```yaml
# Correct
parent:
  child1: value1
  child2: value2

# Incorrect (mismatched indentation)
parent:
  child1: value1
 child2: value2
```

---

## 2. Key-Value Pairs

The most basic element of a YAML file is a key followed by a colon, a space, and a value.

```yaml
label: "mm10"
threads: 8
enable_feature: true
```

---

## 3. Lists (Sequences)

Lists are used when an order matters or when you have multiple values for a single key (like a list of libraries). They start with a dash `-` followed by a space.

```yaml
selected_chromosomes:
  - chr1
  - chr2
  - chrX
```

---

## 4. Dictionaries (Mappings)

Dictionaries are sets of key-value pairs nested under a parent key.

```yaml
genome:
  label: "mm10"
  fasta_path: "refs/genome.fa"
```

---

## 5. Comments

Anything following a `#` is a comment and will be ignored by the pipeline. Use them to document your settings!

```yaml
# This is a comment
defaults:
  strandedness: 0 # 0 means unstranded
```

---

## 6. Common Pitfalls

1. **Missing Space after Colon**: `key:value` is invalid; it must be `key: value`.
2. **Tabs**: If your editor inserts tabs, Snakemake will throw an error. Use a modern editor like VS Code or Vim configured for spaces.
3. **Quotes**: Usually not required unless the value contains special characters or starts with a dash.

---

## 7. Useful Resources

If you want to dive deeper into YAML, check out these excellent resources:

- **[Official YAML Specification](https://yaml.org/spec/1.2.2/)**: The definitive technical guide.
- **[Learn X in Y Minutes (YAML)](https://learnxinyminutes.com/docs/yaml/)**: A great interactive primer.
- **[YAML Validator](https://yamlvalidator.com/)**: Paste your code here to check for syntax errors before running the pipeline.
- **[Wikipedia: YAML](https://en.wikipedia.org/wiki/YAML)**: Good historical and conceptual overview.
