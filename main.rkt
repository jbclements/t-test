#lang typed/racket
(require math/statistics math/special-functions)

(provide student-t-test welch-t-test)

;; given two sequences of numbers, return the two-tailed p-value
;; of the student's t-test. If #:t-statistic is true, return instead
;; the t-statistic. Note that if the two variances are not known
;; to be equal, the welch's t-test (below) is probably a better
;; choice.
(define (student-t-test [S1 : (Sequenceof Real)]
                        [S2 : (Sequenceof Real)]
                        #:t-statistic [t-statistic? : Boolean #f])
  : Real
  (define fun-name 'student-p-value)
  (match-define (list variance1 N1 variance2 N2)
    (compute-variances-and-lengths S1 S2 fun-name))
  (define N1-1 (- N1 1))
  (define N2-1 (- N2 1))
  (define pooled-stdev
    ;; at least one variance positive, both counts positive, so
    ;; numerator should be positive, denominator also positive.
    (sqrt (cast
           (/ (+ (* N1-1 variance1) (* N2-1 variance2))
              (+ N1-1 N2-1))
           Positive-Real)))

  (define t-statistic
    (/ (- (mean S1) (mean S2))
       (* pooled-stdev (sqrt (+ (/ 1 N1) (/ 1 N2))))))

  (cond
    [t-statistic?
     t-statistic]
    [else
     (define degrees-of-freedom
       ;; can't both be zero, by earlier test
       (cast (+ N1-1 N2-1) Positive-Integer))
  
     (t-statistic->p-value t-statistic degrees-of-freedom)]))

;; compute the two-tailed p-value of a Welch's t-test
;; on the two given sequences of real numbers.
(define (welch-t-test [S1 : (Sequenceof Real)]
                      [S2 : (Sequenceof Real)]
                      #:t-statistic [t-statistic? : Boolean #f])
  : Real
  (define fun-name 'student-p-value)
  (match-define (list variance-1 N1 variance-2 N2)
    (compute-variances-and-lengths S1 S2 fun-name))
  ;; N1 & N2 > 0, variances >= 0, therefore these are >= 0.
  ;; Type for / (I claim) doesn't cover this case:
  (define v1/n1 (cast (/ variance-1 N1) Nonnegative-Real))
  (define v2/n2 (cast (/ variance-2 N2) Nonnegative-Real))
 
  (define t-statistic
    (/ (- (mean S1) (mean S2))
       (sqrt (+ v1/n1 v2/n2))))

  ;; the t-statistic is distributed according to the t-distribution
  ;; with a number of degrees of freedom approximated as follows:

  (cond
    [t-statistic? t-statistic]
    [else
     ;; uses the Welch-Satterthwaite approximation:
     (define degrees-of-freedom
       ;; once again, Types for / aren't rich enough (I claim)
       (cast
        (/ (sqr (+ v1/n1 v2/n2))
           (+ (/ (sqr variance-1) (* (sqr N1) (sub1 N1)))
              (/ (sqr variance-2) (* (sqr N2) (sub1 N2)))))
        Positive-Real))


     (t-statistic->p-value t-statistic degrees-of-freedom)]))

(define (compute-variances-and-lengths [S1 : (Sequenceof Real)]
                                       [S2 : (Sequenceof Real)]
                                       [fun-name : Symbol])
  : (List Nonnegative-Real Positive-Integer
          Nonnegative-Real Positive-Integer)
  
  (define N1 (let ([x (sequence-length S1)])
               (cond [(= x 0)
                      (raise-argument-error fun-name
                                            "sequence of length > 0"
                                            0 S1 S2)]
                     [else x])))
  (define N2 (let ([x (sequence-length S2)])
               (cond [(= x 0)
                      (raise-argument-error fun-name
                                            "sequence of length > 0"
                                            1 S1 S2)]
                     [else x])))
  (define variance-1 (variance S1 #:bias #t))
  (define variance-2 (variance S2 #:bias #t))
  (when (and (= variance-1 0)
             (= variance-2 0))
    (raise-argument-error fun-name
                          "at least one sequence with nonzero variance"
                          1 S1 S2))
  (list variance-1 N1 variance-2 N2))

;; given a t-statistic drawn from the t-test distribution and the
;; degrees of freedom of the distribution, return the corresponding
;; two-tailed p-value, also known as 2*(1 - CDF(t)).
(define (t-statistic->p-value [t-statistic : Real]
                              [degrees-of-freedom : Positive-Real])
  : Flonum
  
  ;; x is the parameter to the incomplete beta function
  (define x (/ degrees-of-freedom
               (+ (sqr t-statistic) degrees-of-freedom)))
  ;; two-tailed:
  (beta-inc (/ degrees-of-freedom 2) 1/2 x #f #t)
  )
 
(module+ test
  (require typed/rackunit)

  ;; tests from wikipedia:
  (define wikipedia-l1 (list 30.02 29.99 30.11 29.97 30.01 29.99))
  (define wikipedia-l2 (list 29.89 29.93 29.72 29.98 30.02 29.98))
  
  (check-= (welch-t-test wikipedia-l1 wikipedia-l2) 0.09077 1e-4)
  (check-= (student-t-test wikipedia-l1 wikipedia-l2) 0.07857 1e-4)
  (check-= (student-t-test wikipedia-l1 wikipedia-l2 #:t-statistic #t)
   1.959 1e-3)
  
  ;; tests from rosettacode.com, via Tim Brown
  (check-=
   (welch-t-test
    (list 27.5 21.0 19.0 23.6 17.0 17.9 16.9 20.1 21.9 22.6 23.1
          19.6 19.0 21.7 21.4)
    (list 27.1 22.0 20.8 23.4 23.4 23.5 25.8 22.0 24.8 20.2 21.9
          22.1 22.9 20.5 24.4))
   0.021378001462867
   1e-14)

  ;; check that it's symmetric?
  (check-=
   (welch-t-test
    (list 27.1 22.0 20.8 23.4 23.4 23.5 25.8 22.0 24.8 20.2 21.9
          22.1 22.9 20.5 24.4)
    (list 27.5 21.0 19.0 23.6 17.0 17.9 16.9 20.1 21.9 22.6 23.1
          19.6 19.0 21.7 21.4))
   0.021378001462867
   1e-14)

  (check-=
   (welch-t-test
    (list 17.2 20.9 22.6 18.1 21.7 21.4 23.5 24.2 14.7 21.8)
    (list 21.5 22.8 21.0 23.0 21.6 23.6 22.5 20.7 23.4 21.8
          20.7 21.7 21.5 22.5 23.6 21.5 22.5 23.5 21.5 21.8))
   0.148841696605327
   1e-14)

  (check-=
   (welch-t-test
    (list 19.8 20.4 19.6 17.8 18.5 18.9 18.3 18.9 19.5 22.0)
    (list 28.2 26.6 20.1 23.3 25.2 22.1 17.7 27.6 20.6 13.7
          23.2 17.5 20.6 18.0 23.9 21.6 24.3 20.4 24.0 13.2))
   0.0359722710297968
   1e-14)

  (check-=
   (welch-t-test
    (list 30.02 29.99 30.11 29.97 30.01 29.99)
    (list 29.89 29.93 29.72 29.98 30.02 29.98))
   0.090773324285671
   1e-14)
 
  (check-=
   (welch-t-test (list 3.0 4.0 1.0 2.1)
            (list 490.2 340.0 433.9))
   0.0107515611497845
   1e-14)

  (check-=
   (welch-t-test (list 0.010268 0.000167 0.000167)
            (list 0.159258 0.136278 0.122389))
   0.00339907162713746
   1e-14)

  (check-=
   (welch-t-test (list (/ 1.0 15) (/ 10.0 62.0))
            (list (/ 1.0 10) (/ 2 50.0)))
   0.52726574965384
   1e-14)

  (check-=
   (welch-t-test (list (/ 9 23.0) (/ 21 45.0) (/ 0 38.0))
            (list (/ 0 44.0) (/ 42 94.0) (/ 0 22.0)))
   0.545266866977794
   1e-14)

  (check-exn
   #px"at least one sequence with nonzero variance"
   (λ () (welch-t-test '(3 3 3 3) '(4 4 4 4))))

  (check-exn
   #px"sequence of length > 0"
   (λ () (welch-t-test '() '(4 5 4 4))))

  (check-exn
   #px"sequence of length > 0"
   (λ () (welch-t-test '(4 5 4 4) '())))
)
