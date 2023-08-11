#!/bin/bash

# Convert README.md to index.html using pandoc
pandoc README.md -s --metadata pagetitle="Photo Ambiguities" -o index.html

# Link the container.css to index.html
sed -i 's/<\/head>/<link rel="stylesheet" href="container.css"><\/head>/' index.html

# Wrap the body content inside the .container div
sed -i 's/<body>/<body>\n<div class="container">/' index.html
sed -i 's/<\/body>/<\/div>\n<\/body>/' index.html

echo "Processing complete!"
