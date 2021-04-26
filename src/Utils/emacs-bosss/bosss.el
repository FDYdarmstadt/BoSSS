;;; bosss.el --- Support for the BoSSS solver

;; Copyright (C) 2019 Dario Klingenberg

;; Author: Dario Klingenberg <dario.klingenberg@web.de>
;; Created: 30 Aug 2019
;; Keywords: languages
;; Homepage: http://github.com/dakling/emacs-bosss
;; Version: 0.2

;; This file is not part of GNU Emacs.

;; This file is free software
;; along with this file.  If not, see <https://www.gnu.org/licenses/>.

;;; Code:

;;; repl stuff

;; (defvar bosss-repl-path "/usr/bin/bosss-console")
(defvar bosss-repl-path "bosss-console")
;; (defvar bosss-repl-path "/usr/local/bin/sdb")
;; (defvar bosss-repl-path "/usr/bin/mono")

(defvar bosss-repl-arguments nil)
;; (defvar bosss-repl-arguments (list bosss-pad-path))

(defvar bosss-repl-mode-map (make-sparse-keymap))

(defvar bosss-repl-prompt-regexp "^\\(?:\\[[^@]+@[^@]+\\]\\)")

;; (defvar bosss-repl-prompt-regexp ">")

(defvar bosss--block-beginning-mark  "==============")

(defvar bosss--block-end-mark  "**************")

(defun bosss-repl-load-my-assembly ()
  "load the personal project into the repl"
  (interactive)
  (when bosss-path-reference
    (mapcar (lambda (reference) (comint-send-string "*bosss*" (concat "LoadAssembly(\"" reference "\")\n"))) bosss-path-reference)))

(defun bosss-repl-start-bosss-pad ()
  "only makes sense if bossspad is wrapped inside a debugger"
  (interactive)
  ;; (comint-send-string
  ;;  "*bosss*"
  ;;  "args --simpleconsole \n")
  ;; (comint-send-string
  ;;  "*bosss*"
  ;;  (concat "run " bosss-pad-path "\n"))
  )

(defun run-bosss-repl ()
  "start the repl in a new or existing buffer"
  (interactive)
  (let* ((bosss-program bosss-repl-path)
	 (buffer (comint-check-proc "bosss")))
    (display-buffer
     (get-buffer-create (or buffer "*bosss*")))
    (unless buffer
      (apply 'make-comint-in-buffer "bosss" buffer
       bosss-program bosss-repl-arguments)
      (with-current-buffer (or buffer "*bosss*")
        (bosss-repl-mode))))
  (bosss-repl-start-bosss-pad))

(defun bosss--repl-initialize ()
  (setq comint-process-echoes nil)
  (setq comint-input-sender-no-newline nil)
  (setq comint-use-prompt-regexp t))

;;;###autoload
(define-derived-mode bosss-repl-mode comint-mode "bosss-repl"
  (setq comint-prompt-regexp bosss-repl-prompt-regexp)
  (setq comint-prompt-read-only t)
  (set (make-local-variable 'paragraph-separate) "\\'")
  ;; (set (make-local-variable 'font-lock-defaults) '(bosss-repl-font-lock-keywords t))
  (set (make-local-variable 'paragraph-start) bosss-repl-prompt-regexp))

(add-hook 'bosss-repl-mode-hook 'bosss--repl-initialize)

(defun run-bosss-repl-other-window ()
  "start the repl in another window"
    (interactive)
    (switch-to-buffer-other-window "*bosss*")
    (run-bosss-repl))

;; (defun bosss-repl-send-region (beg end)
;;   (interactive "r")
;;   (comint-send-region "*bosss*" beg end))

(defun bosss--format-input (input-string)
  "formats the input as a single line, which is better handled by comint"
  ;; remove linebreaks
 (replace-regexp-in-string
                       "\n"
                       ""
  ;; remove comments first
                       (replace-regexp-in-string
                        (rx (or (and "\*" (*? anything) "*/") (and "//" (*? anything) eol)))
                        ""
                         input-string)))

(defun bosss-repl-send-string (string-to-send)
  "send the active region to comint"
  (comint-send-string "*bosss*"
                      (concat
                       (bosss--format-input
                        string-to-send)
                       "\n")))

(defun bosss-repl-send-region (beg end)
  "send the active region to comint"
  (interactive "r")
  (bosss-repl-send-string (buffer-substring beg end)))

(defun bosss-repl-send-current-field ()
  "send the current field to comint"
  (interactive)
  (save-excursion
   (search-backward bosss--block-beginning-mark)
   (forward-line 1)
   (let ((beg (point)))
     (search-forward bosss--block-end-mark)
     (move-end-of-line 0)
     (bosss-repl-send-region beg (point)))))

(defun bosss-repl-quit ()
  "send the current field to comint"
  (interactive)
  (bosss-repl-send-string "quit"))

;; not working; we need to check if each command has finished!
;; (defun bosss-repl-send-buffer ()
;;   (interactive)
;;   (save-excursion
;;     (goto-line 0)
;;     (while (search-forward bosss--block-beginning-mark nil t)
;;       (bosss-eval-and-next-field))))


;; worksheet stuff

;;;###autoload
(define-derived-mode bosss-mode csharp-mode "bosss"
  ;; (setq bosss-highlights
  ;; 	'("==============" . font-lock-comment-face))
  ;; ;; (setq bosss-highlights
  ;; ;; 	'(("==============" . font-lock-comment-face)
  ;; ;;     ("**************" . font-lock-comment-face)))
  ;; (setq font-lock-defaults (cons (list bosss-highlights (car csharp-font-lock-keywords)) (cdr font-lock-defaults)))
  )

(defvar bosss-mode-map (make-sparse-keymap))

(defun bosss-create-new-field ()
  (interactive)
  (search-forward bosss--block-end-mark nil t)
  (newline)
  (insert (concat "// " bosss--block-beginning-mark))
  (newline)
  (newline)
  (newline)
  (insert (concat "// " bosss--block-end-mark))
  (forward-line -2))

(defun bosss-next-field ()
  (interactive)
  (if (search-forward bosss--block-beginning-mark nil t)
      (forward-line 1)
    (message "Already on last field")))

(defun bosss-previous-field ()
  (interactive)
  (if (search-backward bosss--block-beginning-mark nil t 2)
      (forward-line 1)
    (message "Already on first field")))

(defun bosss-eval-and-next-field ()
  (interactive)
  (bosss-repl-send-current-field)
  (bosss-next-field))

;;ugly but working
(defun bosss-comment-all-separators ()
  "Comment out all separators"
  (interactive)
   (replace-regexp "^===" "// ===" nil (point-min) (point-max))
   (replace-regexp "^\\*\\*\\*" "// \*\*\*" nil (point-min) (point-max)))

;;ugly but working
(defun bosss-uncomment-all-separators ()
  "Uncomment all separators"
  (interactive)
   (replace-regexp "^// ===" "===" nil (point-min) (point-max))
   (replace-regexp "^// \\*\\*\\*" "\*\*\*" nil (point-min) (point-max)))

(defun bosss-get-most-recent-deploy-directory ()
  (with-current-buffer "*bosss*"
    (progn
     (goto-char (point-max))
     (search-backward "Deployment directory:" nil t)
     (end-of-line)
     (thing-at-point 'filename t))))

;; a bit of a HACK
(defun bosss-get-most-recent-pid ()
  (with-current-buffer "*bosss*"
    (progn
     (goto-char (point-max))
     (search-backward-regexp "^[0-9]" nil t)
     (end-of-line)
     (thing-at-point 'word t))))

;; TODO define text object for a field

(provide 'bosss)

;;; bosss.el ends here
