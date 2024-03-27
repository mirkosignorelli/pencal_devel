# step 1
lmmfit = function(x, ...) UseMethod('lmmfit', x)
mlpmmfit = function(x, ...) UseMethod('mlpmmfit', x)

# step 2
ranefs = function(x, ...) UseMethod('ranefs', x)

# step 3
prclmm = function(x, ...) UseMethod('prclmm', x)
prcmlpmm = function(x, ...) UseMethod('prcmlpmm', x)
sprclmm = function(x, ...) UseMethod('sprclmm', x)
sprcmlpmm = function(x, ...) UseMethod('sprcmlpmm', x)
