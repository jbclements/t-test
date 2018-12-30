#lang info

(define collection "t-test")
(define version "1.0")
(define deps '("base"
               "math-lib"
               "typed-racket-lib"))
(define build-deps '("racket-doc"
                     "rackunit-typed"
                     "scribble-lib"))

(define scribblings '(("t-test.scrbl" ())))