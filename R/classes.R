# step 1
lmmfit = function(x, ...) UseMethod('lmmfit', x)
mlpmmfit = function(x, ...) UseMethod('mlpmmfit', x)

# step 2
lmmsum = function(x, ...) UseMethod('lmmsum', x)
mlpmmsum = function(x, ...) UseMethod('mlpmmsum', x)

# step 3
prclmm = function(x, ...) UseMethod('prclmm', x)
prcmlpmm = function(x, ...) UseMethod('prcmlpmm', x)
sprclmm = function(x, ...) UseMethod('sprclmm', x)
sprcmlpmm = function(x, ...) UseMethod('sprcmlpmm', x)
