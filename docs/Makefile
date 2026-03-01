# docs/Makefile — minimal Sphinx build wrapper
#
# Usage (from the docs/ directory):
#   make html        — build HTML docs into docs/_build/html/
#   make clean       — remove the _build directory
#   make livehtml    — auto-rebuild on file changes (requires sphinx-autobuild)

SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help html clean livehtml

help:
	@echo "make html      — build HTML documentation"
	@echo "make clean     — remove _build/"
	@echo "make livehtml  — live-reload server (requires sphinx-autobuild)"

html:
	$(SPHINXBUILD) -b html $(SOURCEDIR) $(BUILDDIR)/html
	@echo ""
	@echo "Build finished. Open $(BUILDDIR)/html/index.html in your browser."

clean:
	rm -rf $(BUILDDIR)

livehtml:
	sphinx-autobuild $(SOURCEDIR) $(BUILDDIR)/html
