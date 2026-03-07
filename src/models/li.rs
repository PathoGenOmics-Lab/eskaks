use crate::codon::INVALID_CODON;

// --- Constantes del modelo Li ---

/// Codigo genetico: 0 = estandar, != 0 = mitocondrial.
const CODE_MT: i32 = 0;

/// Fraccion 1/3 usada en ponderacion de comparaciones de codones con 2 diferencias.
const ONE_THIRD: f64 = 1.0 / 3.0;

/// Denominador minimo para evitar division por cero en el modelo Li.
const LI_EPSILON: f64 = 1e-15;

/// Matriz de similitud de aminoacidos (Li 1993).
static MAT: [[f64; 19]; 19] = [
    [0.382,0.382,0.343,0.382,0.382,0.382,0.382,0.128,0.040,0.128,0.040,0.128,0.040,0.128,0.040,0.128,0.343,0.128,0.040],
    [0.382,0.382,0.128,0.343,0.343,0.343,0.343,0.128,0.040,0.128,0.040,0.128,0.040,0.128,0.040,0.128,0.128,0.040,0.040],
    [0.343,0.128,0.343,0.382,0.382,0.382,0.343,0.128,0.040,0.128,0.128,0.343,0.128,0.343,0.128,0.343,0.343,0.128,0.040],
    [0.382,0.343,0.382,0.343,0.343,0.343,0.343,0.343,0.040,0.343,0.343,0.382,0.343,0.382,0.343,0.382,0.382,0.382,0.343],
    [0.382,0.343,0.382,0.343,0.382,0.382,0.382,0.343,0.040,0.343,0.128,0.343,0.128,0.128,0.128,0.343,0.343,0.128,0.040],
    [0.382,0.343,0.382,0.343,0.382,0.382,0.382,0.343,0.040,0.343,0.128,0.343,0.128,0.128,0.040,0.128,0.128,0.128,0.040],
    [0.382,0.343,0.343,0.343,0.382,0.382,0.382,0.343,0.040,0.343,0.128,0.343,0.128,0.128,0.128,0.128,0.343,0.128,0.040],
    [0.128,0.128,0.128,0.343,0.343,0.343,0.343,0.343,0.040,0.343,0.128,0.343,0.128,0.343,0.128,0.343,0.343,0.128,0.040],
    [0.040,0.040,0.040,0.040,0.040,0.040,0.040,0.040,0.040,0.382,0.382,0.382,0.343,0.343,0.343,0.128,0.128,0.343,0.128],
    [0.128,0.128,0.128,0.343,0.343,0.343,0.343,0.343,0.382,0.040,0.040,0.128,0.128,0.040,0.128,0.040,0.040,0.040,0.040],
    [0.040,0.040,0.128,0.343,0.128,0.128,0.128,0.128,0.382,0.040,0.343,0.343,0.343,0.343,0.128,0.128,0.128,0.128,0.128],
    [0.128,0.128,0.343,0.382,0.343,0.343,0.343,0.343,0.382,0.128,0.343,0.343,0.343,0.343,0.343,0.128,0.128,0.343,0.343],
    [0.040,0.040,0.128,0.343,0.128,0.128,0.128,0.128,0.343,0.128,0.343,0.343,0.343,0.382,0.343,0.343,0.343,0.343,0.343],
    [0.128,0.128,0.343,0.382,0.128,0.128,0.128,0.343,0.343,0.040,0.343,0.343,0.382,0.343,0.382,0.128,0.128,0.343,0.343],
    [0.040,0.040,0.128,0.343,0.128,0.040,0.128,0.128,0.343,0.128,0.128,0.343,0.343,0.382,0.382,0.343,0.382,0.382,0.343],
    [0.128,0.128,0.343,0.382,0.343,0.128,0.128,0.343,0.128,0.040,0.128,0.128,0.343,0.128,0.343,0.343,0.343,0.382,0.382],
    [0.343,0.128,0.343,0.382,0.343,0.128,0.343,0.343,0.128,0.040,0.128,0.128,0.343,0.128,0.382,0.343,0.382,0.343,0.128],
    [0.128,0.040,0.128,0.382,0.128,0.128,0.128,0.128,0.343,0.040,0.128,0.343,0.343,0.343,0.382,0.382,0.343,0.343,0.343],
    [0.040,0.040,0.040,0.343,0.040,0.040,0.040,0.040,0.128,0.040,0.128,0.343,0.343,0.343,0.343,0.382,0.128,0.343,0.382],
];

// --- Funciones biologicas internas ---

fn codon_to_index(cod: &[char; 3]) -> usize {
    let base_to_num = |c: char| match c { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'U' => 3, _ => 4 };
    let n1 = base_to_num(cod[0]);
    let n2 = base_to_num(cod[1]);
    let n3 = base_to_num(cod[2]);
    if n1 == 4 || n2 == 4 || n3 == 4 { INVALID_CODON as usize } else { 16 * n1 + 4 * n2 + n3 }
}

fn decode_codon(n: usize) -> [char; 3] {
    let b1 = n / 16;
    let r1 = n % 16;
    let b2 = r1 / 4;
    let b3 = r1 % 4;
    let map = |x: usize| match x { 0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T', _ => 'X' };
    [map(b1), map(b2), map(b3)]
}

fn categorize_site(c1: char, c2: char, c3: char, i: i32) -> usize {
    if i == 3 {
        if CODE_MT == 0 {
            if (c1 == 'A' && c2 == 'T' && c3 == 'G') || (c1 == 'T' && c2 == 'G' && (c3 == 'A' || c3 == 'G')) {
                return 0;
            }
        }
        if c2 == 'C' { return 2; }
        if (c1 == 'C' && c2 == 'T') || (c1 == 'G' && c2 == 'T') || (c1 == 'G' && c2 == 'G') || (c1 == 'C' && c2 == 'G') {
            return 2;
        }
        return 1;
    } else if i == 1 {
        if (c1 == 'C' && c2 == 'T' && (c3 == 'A' || c3 == 'G')) || (c1 == 'T' && c2 == 'T' && (c3 == 'A' || c3 == 'G')) {
            return 1;
        }
        if CODE_MT == 0 {
            if (c1 == 'A' && c2 == 'G' && (c3 == 'A' || c3 == 'G')) || (c1 == 'C' && c2 == 'G' && (c3 == 'A' || c3 == 'G')) {
                return 1;
            }
        }
        return 0;
    }
    0
}

fn classify_mutation(nt1: char, nt2: char) -> char {
    if nt1 == nt2 {
        'S'
    } else {
        match (nt1, nt2) {
            ('A','C')|('A','T')|('C','A')|('G','C')|('G','T')|('T','A')|('C','G')|('T','G') => 'v',
            ('A','G')|('C','T')|('G','A')|('T','C') => 'i',
            _ => 'E',
        }
    }
}

fn fill_aa(aa: &mut [usize; 64]) {
    aa[0]=17;aa[1]=16;aa[2]=17;aa[3]=16; aa[4]=13;aa[5]=13;aa[6]=13;aa[7]=13;
    if CODE_MT!=0 {aa[8]=0;} else {aa[8]=18;}
    aa[9]=14;
    if CODE_MT!=0 {aa[10]=0;} else {aa[10]=18;}
    aa[11]=14; if CODE_MT!=0 {aa[12]=5;} else {aa[12]=7;}
    aa[13]=7;aa[14]=5;aa[15]=7; aa[16]=15;aa[17]=4;aa[18]=15;aa[19]=4;aa[20]=9;aa[21]=9;aa[22]=9;aa[23]=9;
    aa[24]=18;aa[25]=18;aa[26]=18;aa[27]=18; aa[28]=6;aa[29]=6;aa[30]=6;aa[31]=6;
    aa[32]=19;aa[33]=20;aa[34]=19;aa[35]=20; aa[36]=11;aa[37]=11;aa[38]=11;aa[39]=11;
    aa[40]=12;aa[41]=12;aa[42]=12;aa[43]=12; aa[44]=8;aa[45]=8;aa[46]=8;aa[47]=8;
    aa[48]=0;aa[49]=3;aa[50]=0;aa[51]=3; aa[52]=14;aa[53]=14;aa[54]=14;aa[55]=14;
    if CODE_MT!=0 {aa[56]=2;} else {aa[56]=0;}
    aa[57]=10;aa[58]=2;aa[59]=10;aa[60]=6;aa[61]=1;aa[62]=6;aa[63]=1;
}

fn special_titv_adjust(ci1:char,ci2:char,ci3:char,cj1:char,cj2:char,cj3:char,ti:&mut[f64;3],tv:&mut[f64;3],poids:f64) {
    if ci1=='C'&&ci2=='G'&&ci3=='A'&&cj1=='T'&&cj2=='G'&&cj3=='A' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; } if ci1=='C'&&ci2=='G'&&ci3=='G'&&cj1=='T'&&cj2=='G'&&cj3=='G' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; }
    if ci1=='A'&&ci2=='G'&&ci3=='G'&&cj1=='G'&&cj2=='G'&&cj3=='G' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; } if ci1=='A'&&ci2=='G'&&ci3=='A'&&cj1=='G'&&cj2=='G'&&cj3=='A' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; }
    if ci1=='T'&&ci2=='G'&&ci3=='A'&&cj1=='C'&&cj2=='G'&&cj3=='A' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; } if ci1=='T'&&ci2=='G'&&ci3=='G'&&cj1=='C'&&cj2=='G'&&cj3=='G' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; }
    if ci1=='G'&&ci2=='G'&&ci3=='G'&&cj1=='A'&&cj2=='G'&&cj3=='G' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; } if ci1=='G'&&ci2=='G'&&ci3=='A'&&cj1=='A'&&cj2=='G'&&cj3=='A' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; }
    if ci1=='C'&&ci2=='G'&&ci3=='A'&&cj1=='A'&&cj2=='G'&&cj3=='A' { tv[1]-=poids;ti[1]+=poids; } if ci1=='A'&&ci2=='G'&&ci3=='A'&&cj1=='C'&&cj2=='G'&&cj3=='A' { tv[1]-=poids;ti[1]+=poids; }
    if ci1=='C'&&ci2=='G'&&ci3=='G'&&cj1=='A'&&cj2=='G'&&cj3=='G' { tv[1]-=poids;ti[1]+=poids; } if ci1=='A'&&ci2=='G'&&ci3=='G'&&cj1=='C'&&cj2=='G'&&cj3=='G' { tv[1]-=poids;ti[1]+=poids; }
}

fn special_titv_adjust_pos3(ci1:char,ci2:char,ci3:char,cj1:char,cj2:char,cj3:char,ti:&mut[f64;3],tv:&mut[f64;3],poids:f64) {
    if ci1=='A'&&ci2=='T'&&ci3=='A'&&cj1=='A'&&cj2=='T'&&cj3=='T' { tv[1]-=poids;ti[1]+=poids; } if ci1=='A'&&ci2=='T'&&ci3=='T'&&cj1=='A'&&cj2=='T'&&cj3=='A' { tv[1]-=poids;ti[1]+=poids; }
    if ci1=='A'&&ci2=='T'&&ci3=='A'&&cj1=='A'&&cj2=='T'&&cj3=='C' { tv[1]-=poids;ti[1]+=poids; } if ci1=='A'&&ci2=='T'&&ci3=='C'&&cj1=='A'&&cj2=='T'&&cj3=='A' { tv[1]-=poids;ti[1]+=poids; }
    if ci1=='A'&&ci2=='T'&&ci3=='A'&&cj1=='A'&&cj2=='T'&&cj3=='G' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; } if ci1=='A'&&ci2=='T'&&ci3=='G'&&cj1=='A'&&cj2=='T'&&cj3=='A' { ti[1]-=0.5*poids;tv[1]+=0.5*poids; }
}

fn count_substitutions_1diff(cod1:&[char;3],cod2:&[char;3],poids:f64,ti:&mut[f64;3],tv:&mut[f64;3],l:&mut[f64;3]) {
    let ci1=cod1[0];let ci2=cod1[1];let ci3=cod1[2];
    let cj1=cod2[0];let cj2=cod2[1];let cj3=cod2[2];
    for i in 0..3 {
        if cod1[i]!=cod2[i] {
            let site=categorize_site(ci1,ci2,ci3,(i+1)as i32); l[site]+=0.5*poids;
            let site2=categorize_site(cj1,cj2,cj3,(i+1)as i32); l[site2]+=0.5*poids;
            let a=cod1[i];let b=cod2[i];
            if classify_mutation(a,b)=='i' { ti[site]+=0.5*poids;ti[site2]+=0.5*poids; } else { tv[site]+=0.5*poids;tv[site2]+=0.5*poids; }
            if CODE_MT==0 { if (ci2=='T'&&cj2=='T')||(ci2=='G'&&cj2=='G') { if i==0 { special_titv_adjust(ci1,ci2,ci3,cj1,cj2,cj3,ti,tv,poids); } if i==2 { special_titv_adjust_pos3(ci1,ci2,ci3,cj1,cj2,cj3,ti,tv,poids); } } }
        }
    }
}

fn count_substitutions_2diff(cod1:&[char;3],cod2:&[char;3],ti:&mut[f64;3],tv:&mut[f64;3],l:&mut[f64;3],aa:&[usize;64],rl:&[Vec<f64>],pos_diff_flags:&[i32;3]) {
    let mut diff_indices = [0; 2]; let mut current_diff_count = 0;
    for k in 0..3 { if pos_diff_flags[k] == 1 { if current_diff_count < 2 { diff_indices[current_diff_count] = k; } current_diff_count += 1; } }
    if current_diff_count != 2 { return; }
    let d1_idx = diff_indices[0]; let d2_idx = diff_indices[1];
    let mut codint1 = *cod1; codint1[d1_idx] = cod2[d1_idx];
    let mut codint2 = *cod1; codint2[d2_idx] = cod2[d2_idx];
    let aa1=aa[codon_to_index(cod1)]; let aa2=aa[codon_to_index(cod2)];
    let aaint1_path1=aa[codon_to_index(&codint1)]; let aaint1_path2=aa[codon_to_index(&codint2)];
    let l_path1 = rl[aa1][aaint1_path1] * rl[aaint1_path1][aa2];
    let l_path2 = rl[aa1][aaint1_path2] * rl[aaint1_path2][aa2];
    let p1 = if (l_path1 + l_path2) != 0.0 { l_path1 / (l_path1 + l_path2) } else { 0.5 };
    let p2 = 1.0 - p1;
    let mut non_diff_pos_plus_1 = 0;
    for k_idx in 0..3 { if pos_diff_flags[k_idx] == 0 { non_diff_pos_plus_1 = k_idx + 1; break; } }
    if non_diff_pos_plus_1 > 0 {
        l[categorize_site(cod1[0],cod1[1],cod1[2],non_diff_pos_plus_1 as i32)]+=ONE_THIRD;
        l[categorize_site(cod2[0],cod2[1],cod2[2],non_diff_pos_plus_1 as i32)]+=ONE_THIRD;
        l[categorize_site(codint1[0],codint1[1],codint1[2],non_diff_pos_plus_1 as i32)]+=ONE_THIRD*p1;
        l[categorize_site(codint2[0],codint2[1],codint2[2],non_diff_pos_plus_1 as i32)]+=ONE_THIRD*p2;
    }
    count_substitutions_1diff(cod1,&codint1,p1,ti,tv,l);
    count_substitutions_1diff(&codint1,cod2,p1,ti,tv,l);
    count_substitutions_1diff(cod1,&codint2,p2,ti,tv,l);
    count_substitutions_1diff(&codint2,cod2,p2,ti,tv,l);
}

fn count_substitutions_3diff(cod1:&[char;3],cod2:&[char;3],ti:&mut[f64;3],tv:&mut[f64;3],l:&mut[f64;3],aa:&[usize;64],rl:&[Vec<f64>]) {
    let mut codint1_paths = [['A'; 3]; 6]; let mut codint2_paths = [['A'; 3]; 6];
    let mut like = [0.0; 6]; let mut path_idx = 0;
    for i in 0..3 { for j in 0..3 {
        if j == i { continue; }
        let mut c2_intermediate = *cod1;
        c2_intermediate[i] = cod2[i]; codint1_paths[path_idx] = c2_intermediate;
        let mut c3_intermediate = c2_intermediate;
        c3_intermediate[j] = cod2[j]; codint2_paths[path_idx] = c3_intermediate;
        let aa_orig = aa[codon_to_index(cod1)];
        let aa_c2_int = aa[codon_to_index(&c2_intermediate)];
        let aa_c3_int = aa[codon_to_index(&c3_intermediate)];
        let aa_final = aa[codon_to_index(cod2)];
        like[path_idx] = rl[aa_orig][aa_c2_int] * rl[aa_c2_int][aa_c3_int] * rl[aa_c3_int][aa_final];
        path_idx += 1;
    } }
    let somli: f64 = like.iter().sum();
    let mut p = [0.0; 6];
    if somli > 0.0 { for i in 0..6 { p[i] = like[i] / somli; } } else { for i in 0..6 { p[i] = 1.0 / 6.0; } }
    for i in 0..6 {
        if p[i] > 0.0 {
            count_substitutions_1diff(cod1, &codint1_paths[i], p[i], ti, tv, l);
            count_substitutions_1diff(&codint1_paths[i], &codint2_paths[i], p[i], ti, tv, l);
            count_substitutions_1diff(&codint2_paths[i], cod2, p[i], ti, tv, l);
        }
    }
}

// --- LiTables: tablas precomputadas con arrays planos ---

/// Tablas de lookup precomputadas para el modelo Li (1993).
/// Usa arrays planos [f64; 4096] (64x64) en lugar de Vec<Vec<f64>>
/// para eliminar la doble indirection de punteros y mejorar la localidad de cache.
pub struct LiTables {
    tl: [[f64; 4096]; 3],
    tti: [[f64; 4096]; 3],
    ttv: [[f64; 4096]; 3],
}

#[inline(always)]
fn idx(i: usize, j: usize) -> usize {
    i * 64 + j
}

impl LiTables {
    /// Construye las tablas de lookup. Incluye la construccion de la matriz rl_li.
    pub fn new() -> Box<Self> {
        // Construir matriz de similitud rl_li
        let mut rl_li = vec![vec![0.0; 21]; 21];
        for i in 2..=20 {
            for j in 1..i {
                rl_li[i][j] = MAT[j - 1][i - 2];
            }
        }
        for i in 1..=20 {
            rl_li[i][i] = 1.0;
            for j in i + 1..=20 {
                rl_li[i][j] = rl_li[j][i];
            }
        }
        let mut minrl = rl_li.get(1).and_then(|row| row.get(1)).copied().unwrap_or(0.01);
        for i in 1..=20 {
            for j in i + 1..=20 {
                if rl_li[i][j] < minrl { minrl = rl_li[i][j]; }
            }
        }
        for i in 0..=20 {
            rl_li[0][i] = minrl;
            rl_li[i][0] = minrl;
        }

        // Allocar tablas con ceros
        let mut tables = Box::new(LiTables {
            tl: [[0.0; 4096]; 3],
            tti: [[0.0; 4096]; 3],
            ttv: [[0.0; 4096]; 3],
        });

        // Llenar tablas (equivalente a prefastlwl)
        let mut aa_li_map = [0usize; 64];
        fill_aa(&mut aa_li_map);

        for i in 0..64 {
            for j in i..64 {
                let mut l_accum = [0.0; 3];
                let mut ti_accum = [0.0; 3];
                let mut tv_accum = [0.0; 3];
                let cod1_chars = decode_codon(i);
                let cod2_chars = decode_codon(j);
                let mut pos_diff_flags = [0i32; 3];
                let mut nbdiff = 0;
                for p in 0..3 {
                    if cod1_chars[p] != cod2_chars[p] {
                        nbdiff += 1;
                        pos_diff_flags[p] = 1;
                    }
                }

                if nbdiff == 0 {
                    for p_idx in 0..3 {
                        let site_type = categorize_site(cod1_chars[0], cod1_chars[1], cod1_chars[2], (p_idx + 1) as i32);
                        l_accum[site_type] += 1.0;
                    }
                } else if nbdiff == 1 {
                    for p_idx in 0..3 {
                        if pos_diff_flags[p_idx] == 0 {
                            let st1 = categorize_site(cod1_chars[0], cod1_chars[1], cod1_chars[2], (p_idx + 1) as i32);
                            l_accum[st1] += 0.5;
                            let st2 = categorize_site(cod2_chars[0], cod2_chars[1], cod2_chars[2], (p_idx + 1) as i32);
                            l_accum[st2] += 0.5;
                        }
                    }
                    count_substitutions_1diff(&cod1_chars, &cod2_chars, 1.0, &mut ti_accum, &mut tv_accum, &mut l_accum);
                } else if nbdiff == 2 {
                    for p_idx in 0..3 {
                        if pos_diff_flags[p_idx] == 0 {
                            let st1 = categorize_site(cod1_chars[0], cod1_chars[1], cod1_chars[2], (p_idx + 1) as i32);
                            l_accum[st1] += 0.5;
                            let st2 = categorize_site(cod2_chars[0], cod2_chars[1], cod2_chars[2], (p_idx + 1) as i32);
                            l_accum[st2] += 0.5;
                            break;
                        }
                    }
                    count_substitutions_2diff(&cod1_chars, &cod2_chars, &mut ti_accum, &mut tv_accum, &mut l_accum, &aa_li_map, &rl_li, &pos_diff_flags);
                } else if nbdiff == 3 {
                    count_substitutions_3diff(&cod1_chars, &cod2_chars, &mut ti_accum, &mut tv_accum, &mut l_accum, &aa_li_map, &rl_li);
                }

                // Llenar simetricamente usando indexacion plana
                for k in 0..3 {
                    tables.tl[k][idx(i, j)] = l_accum[k];
                    tables.tl[k][idx(j, i)] = l_accum[k];
                    tables.tti[k][idx(i, j)] = ti_accum[k];
                    tables.tti[k][idx(j, i)] = ti_accum[k];
                    tables.ttv[k][idx(i, j)] = tv_accum[k];
                    tables.ttv[k][idx(j, i)] = tv_accum[k];
                }
            }
        }

        tables
    }

    /// Calcula Ka y Ks para un par de secuencias codificadas como indices de codones.
    pub fn compute_pair(&self, codon_indices1: &[u8], codon_indices2: &[u8]) -> (f64, f64) {
        let nb_codon = codon_indices1.len();
        let mut l_sum = [0.0; 3];
        let mut ti_sum = [0.0; 3];
        let mut tv_sum = [0.0; 3];

        for k in 0..nb_codon {
            let num1 = codon_indices1[k] as usize;
            let num2 = codon_indices2[k] as usize;
            if num1 == INVALID_CODON as usize || num2 == INVALID_CODON as usize { continue; }
            let flat = idx(num1, num2);
            for t in 0..3 {
                l_sum[t] += self.tl[t][flat];
                ti_sum[t] += self.tti[t][flat];
                tv_sum[t] += self.ttv[t][flat];
            }
        }

        let mut p_k = [0.0; 3];
        let mut q_k = [0.0; 3];
        for k_type in 0..3 {
            if l_sum[k_type] == 0.0 {
                p_k[k_type] = 0.0;
                q_k[k_type] = 0.0;
            } else {
                p_k[k_type] = ti_sum[k_type] / l_sum[k_type];
                q_k[k_type] = tv_sum[k_type] / l_sum[k_type];
            }
        }

        let mut a_k_val = [0.0; 3];
        let mut b_k_val = [0.0; 3];
        for k_type in 0..3 {
            if l_sum[k_type] == 0.0 {
                a_k_val[k_type] = f64::NAN;
                b_k_val[k_type] = f64::NAN;
                continue;
            }
            let denom_a = 1.0 - 2.0 * p_k[k_type] - q_k[k_type];
            let denom_b = 1.0 - 2.0 * q_k[k_type];
            if denom_a <= LI_EPSILON || denom_b <= LI_EPSILON {
                a_k_val[k_type] = f64::NAN;
                b_k_val[k_type] = f64::NAN;
            } else {
                a_k_val[k_type] = -0.5 * denom_a.ln() - 0.25 * denom_b.ln();
                b_k_val[k_type] = -0.25 * denom_b.ln();
            }
        }

        let mut ks_val = f64::NAN;
        let mut ka_val = f64::NAN;
        let l0_s = l_sum[0]; let l2_s = l_sum[1]; let l4_s = l_sum[2];
        let a0_li = a_k_val[0]; let a2_li = a_k_val[1]; let a4_li = a_k_val[2];
        let b0_li = b_k_val[0]; let b2_li = b_k_val[1]; let b4_li = b_k_val[2];

        // Calculo robusto de Ka
        if l0_s > 0.0 && !a0_li.is_nan() && !b0_li.is_nan() {
            if l2_s > 0.0 && !b2_li.is_nan() {
                ka_val = a0_li + (l0_s * b0_li + l2_s * b2_li) / (l0_s + l2_s);
            } else {
                ka_val = a0_li + b0_li;
            }
        }

        // Calculo robusto de Ks
        if (l2_s + l4_s) > 0.0 && !b4_li.is_nan() {
            let num_l2 = if l2_s > 0.0 && !a2_li.is_nan() { l2_s * a2_li } else { 0.0 };
            let den_l2 = if l2_s > 0.0 && !a2_li.is_nan() { l2_s } else { 0.0 };
            let num_l4 = if l4_s > 0.0 && !a4_li.is_nan() { l4_s * a4_li } else { 0.0 };
            let den_l4 = if l4_s > 0.0 && !a4_li.is_nan() { l4_s } else { 0.0 };
            let num = num_l2 + num_l4;
            let den = den_l2 + den_l4;
            if den > 0.0 {
                ks_val = num / den + b4_li;
            }
        }

        if ka_val.is_finite() && ka_val < 0.0 { ka_val = 0.0; }
        if ks_val.is_finite() && ks_val < 0.0 { ks_val = 0.0; }

        (ka_val, ks_val)
    }
}
