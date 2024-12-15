# Variables
POETRY := poetry run
MISS_HIT_STYLE := $(POETRY) mh_style
MISS_HIT_LINT := $(POETRY) mh_lint
MISS_HIT_METRIC := $(POETRY) mh_metric

# Directories to analyze
DIRECTORIES := ./../docs/ ./../test/
# TODO: Add docs examples and src
CONFIG_FILE := ./miss_hit.cfg

# Default target (runs all checks)
.PHONY: all
all: check

# Install project dependencies using Poetry
.PHONY: install
install:
	@echo "Installing dependencies..."
	poetry install

# Format code using mh_style with --fix to apply changes
.PHONY: format
format:
	@echo "Formatting code in $(DIRECTORIES)..."
	$(MISS_HIT_STYLE) $(DIRECTORIES) --fix

# Check code formatting without applying changes
.PHONY: format_check
format_check:
	@echo "Checking code formatting in $(DIRECTORIES)..."
	$(MISS_HIT_STYLE) $(DIRECTORIES)

# Run lint checks using mh_lint
.PHONY: lint
lint:
	@echo "Running lint checks on $(DIRECTORIES)..."
	$(MISS_HIT_LINT) $(DIRECTORIES)

# Run code metrics check with CI mode
.PHONY: metric
metric:
	@echo "Running code metrics for $(DIRECTORIES)..."
	$(MISS_HIT_METRIC) $(DIRECTORIES) --ci

# Run all checks: format_check, lint, and metric
.PHONY: check
check: format lint metric
	@echo "All checks passed successfully."

# Help target to display available commands
.PHONY: help
help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Available targets:"
	@echo "  install         Install project dependencies using Poetry"
	@echo "  format          Format code in directories using mh_style"
	@echo "  format_check    Check code formatting without making changes"
	@echo "  lint            Run lint checks on the specified directories"
	@echo "  metric          Run code metrics in CI mode"
	@echo "  check           Run all checks (format_check, lint, metric)"
	@echo "  clean           Clean up temporary files (if applicable)"
	@echo "  help            Display this help message"
