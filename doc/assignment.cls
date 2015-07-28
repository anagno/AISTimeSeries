%% ============================================================
\NeedsTeXFormat{LaTeX2e}[2005/12/01]
\ProvidesClass{assignment}[2011/09/01 v0.1 DIT.UoA Assignment LaTeX
  document class]

\LoadClass[a4paper,titlepage,twoside,10pt,leqno]{article}[2005/09/16]

%----------------------------------------------------------------------------------------
%	GENERAL PACKAGES
%----------------------------------------------------------------------------------------

% Using include commands without \clearpage
% http://stackoverflow.com/questions/1209065/to-have-no-pagebreak-after-include-in-latex
\RequirePackage{newclude}

\RequirePackage[]{fontspec}
\RequirePackage{xltxtra}

% Τα πακέτα του xelatex που ορίζουν την γλώσσα
\RequirePackage{polyglossia}
\setmainlanguage{english}
\setotherlanguage[variant=mono]{greek}
\setotherlanguage[spelling=new,latesthyphen=true,babelshorthands=true]{german}

% Για υφενώσεις στα ελληνικά
\newcommand{\el}[1]{\begin{greek}#1\end{greek}}

% Για υφενώσεις στα γερμανικά
\newcommand{\de}[1]{\begin{german}#1\end{german}}

% http://tex.stackexchange.com/questions/45598/putting-quotation-marks-in-tex-document
\setmainfont[Ligatures=TeX]{Ubuntu} % Για τις γραμματοσειρές 
\setsansfont[]{Ubuntu} 
\setmonofont[]{Ubuntu} 

% Για τα μαθηματικά μέσω στο κείμενο (π.χ. \eqref)
\RequirePackage{amsmath}
\RequirePackage{amsfonts}

% Landscape pages
\RequirePackage{pdflscape}

% Longtables package
\RequirePackage{longtable}
\RequirePackage{afterpage}

% Για το περιβάλλον inparaenum που είναι οι λίστες μέσα στο κείμενο.
\RequirePackage{paralist} 

% http://stackoverflow.com/questions/5690679/add-appendix-before-a-in-thesis-toc
\RequirePackage[titletoc]{appendix}

% Find the last value of a counter
\RequirePackage{totcount} 
% total number of section
\regtotcounter{section}

%Να δω για να αλλάξω διάφορα
%http://www.mechpedia.gr/wiki/Hyphenation_-_%CE%A5%CF%86%CE%B5%CE%BD%CF%8E%CF%83%CE%B5%CE%B9%CF%82
%http://www.mechpedia.gr/wiki/%CE%95%CE%BB%CE%BB%CE%B7%CE%BD%CE%B9%CE%BA%CE%AC

\RequirePackage{hyperref} % In-text references
% \url command in bibliography

% Για να μην φαίνονται γραμμές γύρω από τις υπερσυνδέσεις μέσα στο κείμενο
\hypersetup{%
pdfborder = {0 0 0},
}

%----------------------------------------------------------------------------------------
%	END GENERAL PACKAGES
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
%	GENERAL ENVIROMENTS
%----------------------------------------------------------------------------------------

\newenvironment{dedication} % Για την αφιέρωση
{
   \cleardoublepage
   \thispagestyle{empty}
   \vspace*{\stretch{1}}
   \hfill\begin{minipage}[t]{0.66\textwidth}
   \raggedleft
}%
{
   \end{minipage}
   \vspace*{\stretch{3}}
   \clearpage
}

\newenvironment{acknowledgements} % Για τις ευχαριστίες
    {\cleardoublepage\thispagestyle{empty}\null\vfill\begin{center}%
    \bfseries \huge{Acknowledgements} \vspace*{10pt} \end{center}}%
    {\vfill\null}

%----------------------------------------------------------------------------------------
%	END GENERAL ENVIROMENTS
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
%	MARGIN OF TEXT PAGES
%----------------------------------------------------------------------------------------

% Τα περιθώρια σελίδας να είναι 3cm και στις τέσσερεις πλευρές (πάνω, κάτω, αριστερά, δεξιά).
\RequirePackage[%
  a4paper,
  top=3cm,   bottom=3cm,     left=3cm,    right=3cm, bindingoffset=0.5cm, 
  head=14pt, headsep=0.50cm, foot=0.75cm, footskip=0.9cm,
  headsep=0.7cm ]{geometry}

%----------------------------------------------------------------------------------------
%	END MARGIN OF TEXT PAGE
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
% Format of Titles for Chapters, Sections, Subsections and Subsubsections.
%----------------------------------------------------------------------------------------

\RequirePackage{titlesec}

%----------------------------------------------------------------------------------------
% End of Format of Titles for Chapters, Sections, Subsections and Subsubsections.
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
% Paragraph and Line Spacing
%----------------------------------------------------------------------------------------

% Single space between lines
\RequirePackage[singlespacing]{setspace} % Lines

% Vertical space after a paragraph should be 3 pt.
\setlength{\parskip}{3pt}

%----------------------------------------------------------------------------------------
% End of Paragraph and Line Spacing
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
% Format of Table of Contents, List of Figures and List of Tables.
%----------------------------------------------------------------------------------------

% Adjust the tilte margins
% http://tex.stackexchange.com/questions/13690/extra-space-for-numbered-titles-in-table-of-contents
\makeatletter
\renewcommand{\l@section}{\@dottedtocline{1}{1.5em}{2.6em}}
\renewcommand{\l@subsection}{\@dottedtocline{2}{4.0em}{3.6em}}
\renewcommand{\l@subsubsection}{\@dottedtocline{3}{7.4em}{4.5em}}
\makeatother

%----------------------------------------------------------------------------------------
% END Format of Table of Contents, List of Figures and List of Tables.
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
%	FIGURES, SUBFIGURES AND CAPTIONS
%----------------------------------------------------------------------------------------

\RequirePackage{tabularx} % Special tables. Επεκτείνει το περιβάλλον tabular
\RequirePackage{multirow} % Για την συγχώνευση πολλαπλών γραμμών στους πίνακες


% Για την αρίθμηση των πινάκων και των σχημάτων ανάλογα με τις ενότητες
\RequirePackage{chngcntr}
\counterwithin{figure}{section}
\counterwithin{table}{section}

\RequirePackage{subcaption}          % Subfigures
% Να κοιτάξω να το φτιάξω σωστά

% Captions
\RequirePackage[%
  justification=centering,
  labelsep=colon,
  figureposition=bottom,
  tableposition=top]{caption}          
% caption fonts should be boldface 10pt.
\renewcommand{\captionfont}{\fontsize{10pt}{12pt}\selectfont\bfseries}

%----------------------------------------------------------------------------------------
%	END FIGURES, SUBFIGURES AND CAPTIONS
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
%	MINITOC
%----------------------------------------------------------------------------------------

%Για την δημιουργία mini περιεχομένων στην αρχή του κάθε κεφαλαίου 

%\RequirePackage{minitoc}
%\dominitoc
%\dominilof
%\dominilot

%\renewcommand{\mtctitle}{}  % Για να μην εμφανίζεται η λέξη Contects
%\mtcsetrules{minitoc}{off}  % Για να εξαφανιστούν κάτι γραμμές
%\renewcommand{\afterminitoc}{\newpage} 

%----------------------------------------------------------------------------------------
%	END MINITOC
%----------------------------------------------------------------------------------------


%----------------------------------------------------------------------------------------
%	BIBLIOGRAPHY PACKAGES
%----------------------------------------------------------------------------------------

% http://tex.stackexchange.com/questions/229638/package-biblatex-warning-babel-polyglossia-detected-but-csquotes-missing
\RequirePackage{csquotes}
\RequirePackage[citestyle=numeric,
				natbib=true,
				maxcitenames=3,
				maxbibnames=100,
				backend=bibtex]{biblatex}
\bibliography{../references}
% http://tex.stackexchange.com/questions/117179/another-problem-with-url-hyphenation-in-bibliography
\setcounter{biburllcpenalty}{9000}

%----------------------------------------------------------------------------------------
%	END BIBLIOGRAPHY PACKAGES
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
%	CODE INCLUSION CONFIGURATION
%----------------------------------------------------------------------------------------

\RequirePackage[outputdir=bin]{minted}
\usemintedstyle{borland}
%Να προσπαθήσω να ορίσω ως προεπιλογή το breaklines=true


%----------------------------------------------------------------------------------------
%	END CODE INCLUSION CONFIGURATION
%----------------------------------------------------------------------------------------

%----------------------------------------------------------------------------------------
%	ALGORITH INCLUSION CONFIGURATION
%----------------------------------------------------------------------------------------

% http://tex.stackexchange.com/questions/124902/algorithm-with-chapter-number
\RequirePackage[linesnumbered,ruled,vlined,resetcount,algosection]{algorithm2e}

% Να τα κοιτάξω καλύτερα αυτά !!!

%----------------------------------------------------------------------------------------
%	ALGORITH CODE INCLUSION CONFIGURATION
%----------------------------------------------------------------------------------------
