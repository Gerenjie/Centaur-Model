import sys

master_doc = 'index'
sys.path.insert(0, '.')

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'numpydoc'
]

autodoc_member_order = 'bysource'
autodoc_default_flags = ['members', 'undoc-members']
autoclass_content = 'both'
