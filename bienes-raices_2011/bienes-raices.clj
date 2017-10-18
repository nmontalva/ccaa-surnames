(require '[clojure.string :as s]
         '[clojure.contrib.combinatorics :as c])

(defn get-rows [s]
  (let [[sep & ss] (.split s "\n")]
    (take-nth 2 (partition-by #{sep} ss))))

;; El archivo pdf primero debe ser convertido a HTML usando pdftohtml
;; y luego hay que sacarle la parte HTML usando:
;; sed -ie 's/<br>//g' archivo.html
;; sed -ie 's/&nbsp;/ /g' archivo.html
;; vim archivo.html
(defn convert [input-filename output-filename]
  (binding [*out* (java.io.FileWriter. output-filename)]
    (println "Numero comunero;Nombres;Apellido paterno;Apellido materno")
    (doseq [[num-comunero nombres apellidos] (get-rows (slurp input-filename))
            :let [apellidos (when apellidos
                              (s/replace apellidos #"\ +" " "))]]
      (println (apply str (interpose ";" (concat [num-comunero nombres]
                                                 (when apellidos
                                                   (.split apellidos " ")))))))))

(defn prob-n [pred input-filename n]
  (let [apellidos (for [[_ _ apellidos] (get-rows (slurp input-filename))
                        :when apellidos]
                    (.split (s/replace apellidos #"\ +" " ") " "))
        conteo (for [as (c/combinations apellidos n)]
                 (if (apply pred as) 1 0))]
    (/ (apply + conteo) (count conteo))))

(defn non-empty-= [& xs]
  (if (some empty? xs)
    false
    (apply = xs)))

(defn ningun-apellido-en-comun [[p1 m1] & as]
  (not-any? #(or (non-empty-= p1 (first %)) (non-empty-= p1 (second %))
                 (non-empty-= m1 (first %)) (non-empty-= m1 (second %))) as))

(defn un-apellido-en-comun [[p1 m1] & as]
  ;; los que tienen dos apellidos en comun deberian contar aqui?
  (every? #(or (non-empty-= p1 (first %)) (non-empty-= p1 (second %))
               (non-empty-= m1 (first %)) (non-empty-= m1 (second %))) as))

(comment ;; otra manera, mas concisa, de escribir lo mismo
(defn un-apellido-en-comun [a1 & as]
  (every? (partial some (set a1)) as)))

(defn dos-apellidos-en-comun-mismo-orden [[p1 m1] & as]
  (every? #(and (non-empty-= p1 (first %)) (non-empty-= m1 (second %))) as))

(defn dos-apellidos-en-comun-sin-orden [[p1 m1] & as]
  (every? #(or (and (non-empty-= p1 (first %)) (non-empty-= m1 (second %)))
               (and (non-empty-= m1 (first %)) (non-empty-= p1 (second %)))) as))

(defn write-prob-table [input-filenames output-filename n]
  (binding [*out* (java.io.FileWriter. output-filename)]
    (println "Comunidad;Prob ningun apellido en comun;Prob un apellido en comun;Prob dos apellidos en comun sin orden;Prob dos apellidos en comun mismo orden")
    (doseq [i input-filenames
            :let [name (s/capitalize (first (.split (slurp i) "\n")))]]
      (println (str name ";"
                    (float (prob-n ningun-apellido-en-comun i n)) ";"
                    (float (prob-n un-apellido-en-comun i n)) ";"
                    (float (prob-n dos-apellidos-en-comun-sin-orden i n)) ";"
                    (float (prob-n dos-apellidos-en-comun-mismo-orden i n)))))))

(defn file-list [start end error-filename]
  (let [erroneous (set (map #(Integer/parseInt %) (.split (slurp error-filename) "\n")))]
    (for [i (range start end)
          :when (not (erroneous i))]
      (str i ".txt"))))

