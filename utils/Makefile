POETRY := poetry run
MISS_HIT_STYLE := $(POETRY) mh_style
MISS_HIT_LINT := $(POETRY) mh_lint
MISS_HIT_METRIC := $(POETRY) mh_metric

# Directories to analyze
DIRECTORIES := ./../docs/ ./../test/ ./../examples

.PHONY: all
all: check

.PHONY: install
install:
	@echo "Installing dependencies..."
	poetry install

.PHONY: format
format:
	@echo "Formatting code in $(DIRECTORIES)..."
	$(MISS_HIT_STYLE) $(DIRECTORIES) --fix

.PHONY: format_check
format_check:
	@echo "Checking code formatting in $(DIRECTORIES)..."
	@FORMAT_OUTPUT=$$($(MISS_HIT_STYLE) $(DIRECTORIES)); \
	echo "$$FORMAT_OUTPUT"; \
	echo "$$FORMAT_OUTPUT" | grep -q ", everything seems fine" && echo "✅ Formatting check passed!" || { echo "❌ Formatting issues found!"; exit 1; }

.PHONY: lint
lint:
	@echo "Running lint checks on $(DIRECTORIES)..."
	@LINT_OUTPUT=$$($(MISS_HIT_LINT) $(DIRECTORIES)); \
	echo "$$LINT_OUTPUT"; \
	echo "$$LINT_OUTPUT" | grep -q ", everything seems fine" && echo "✅ Linting check passed!" || { echo "❌ Linting issues found!"; exit 1; }

.PHONY: metric
metric:
	@echo "Running code metrics for $(DIRECTORIES)..."
	@METRIC_OUTPUT=$$($(MISS_HIT_METRIC) $(DIRECTORIES) --ci); \
	echo "$$METRIC_OUTPUT"; \
	echo "$$METRIC_OUTPUT" | grep -q ", everything seems fine" && echo "✅ Code metrics check passed!" || { echo "❌ Code metrics check failed!"; exit 1; }

.PHONY: check
check: format_check lint metric
	@echo "✅ All checks passed successfully."

.PHONY: ci_checks
ci_checks: install check
	@echo "🎯 CI checks completed successfully."
