import yaml
import jsonschema

with open('config/config.yaml') as f:
    config = yaml.safe_load(f)

with open('workflow/schemas/config.schema.yaml') as f:
    schema = yaml.safe_load(f)

try:
    jsonschema.validate(config, schema)
    print("VALIDATION SUCCESSFUL")
except jsonschema.exceptions.ValidationError as e:
    print("VALIDATION FAILED:", e.message)
