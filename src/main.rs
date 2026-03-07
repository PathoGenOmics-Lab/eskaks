use clap::{Parser, ValueEnum};
use crossbeam::channel::unbounded;
use indicatif::{ProgressBar, ProgressStyle};
use log::info;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;
use std::thread;

// --- CONSTANTES GLOBALES ---
const CODE_MT: i32 = 0;
const AA_ARRAY: [char; 64] = [
    'F', 'F', 'L', 'L', 'L', 'L', 'L', 'L', 'I', 'I', 'I', 'M', 'V', 'V', 'V', 'V',
    'S', 'S', 'S', 'S', 'P', 'P', 'P', 'P', 'T', 'T', 'T', 'T', 'A', 'A', 'A', 'A',
    'Y', 'Y', '*', '*', 'H', 'H', 'Q', 'Q', 'N', 'N', 'K', 'K', 'D', 'D', 'E', 'E',
    'C', 'C', '*', 'W', 'R', 'R', 'R', 'R', 'S', 'S', 'R', 'R', 'G', 'G', 'G', 'G',
];
const SYN_SITE_ARRAY: [usize; 64] = [
    1, 1, 2, 2, 3, 3, 4, 4, 2, 2, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    3, 3, 4, 4, 1, 1, 2, 2, 3, 3, 3, 3,
];
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

// --- FUNCIONES BIOLÓGICAS (SIN CAMBIOS EN SU LÓGICA INTERNA) ---
fn num(cod:&[char;3])->usize {
    let base_to_num=|c:char| match c {'A'=>0,'C'=>1,'G'=>2,'T'=>3,'U'=>3,_=>4};
    let n1=base_to_num(cod[0]); let n2=base_to_num(cod[1]); let n3=base_to_num(cod[2]);
    if n1==4||n2==4||n3==4 {64} else {16*n1+4*n2+n3}
}
fn decode_codon(n:usize)->[char;3] {
    let b1=n/16; let r1=n%16; let b2=r1/4; let b3=r1%4;
    let map=|x:usize| match x {0=>'A',1=>'C',2=>'G',3=>'T',_=>'X'};
    [map(b1),map(b2),map(b3)]
}
fn catsite(c1:char,c2:char,c3:char,i:i32)->usize {
    if i==3 {
        if CODE_MT==0 { if (c1=='A'&&c2=='T'&&c3=='G')||(c1=='T'&&c2=='G'&&(c3=='A'||c3=='G')) { return 0; } }
        if c2=='C' {return 2;}
        if (c1=='C'&&c2=='T')||(c1=='G'&&c2=='T')||(c1=='G'&&c2=='G')||(c1=='C'&&c2=='G') { return 2; }
        return 1;
    } else if i==1 {
        if (c1=='C'&&c2=='T'&&(c3=='A'||c3=='G'))||(c1=='T'&&c2=='T'&&(c3=='A'||c3=='G')) { return 1; }
        if CODE_MT==0 { if (c1=='A'&&c2=='G'&&(c3=='A'||c3=='G'))||(c1=='C'&&c2=='G'&&(c3=='A'||c3=='G')) { return 1; } }
        return 0;
    }
    0
}
fn transf(nt1:char,nt2:char)->char {
    if nt1==nt2 {'S'} else { match (nt1,nt2) { ('A','C')|('A','T')|('C','A')|('G','C')|('G','T')|('T','A')|('C','G')|('T','G')=>'v', ('A','G')|('C','T')|('G','A')|('T','C')=>'i', _=>'E', } }
}
fn fill_aa(aa:&mut[usize;64]) {
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
fn titv1(cod1:&[char;3],cod2:&[char;3],poids:f64,ti:&mut[f64;3],tv:&mut[f64;3],l:&mut[f64;3]) {
    let ci1=cod1[0];let ci2=cod1[1];let ci3=cod1[2]; let cj1=cod2[0];let cj2=cod2[1];let cj3=cod2[2];
    for i in 0..3 {
        if cod1[i]!=cod2[i] {
            let site=catsite(ci1,ci2,ci3,(i+1)as i32); l[site]+=0.5*poids;
            let site2=catsite(cj1,cj2,cj3,(i+1)as i32); l[site2]+=0.5*poids;
            let a=cod1[i];let b=cod2[i];
            if transf(a,b)=='i' { ti[site]+=0.5*poids;ti[site2]+=0.5*poids; } else { tv[site]+=0.5*poids;tv[site2]+=0.5*poids; }
            if CODE_MT==0 { if (ci2=='T'&&cj2=='T')||(ci2=='G'&&cj2=='G') { if i==0 { special_titv_adjust(ci1,ci2,ci3,cj1,cj2,cj3,ti,tv,poids); } if i==2 { special_titv_adjust_pos3(ci1,ci2,ci3,cj1,cj2,cj3,ti,tv,poids); } } }
        }
    }
}
fn titv2(cod1:&[char;3],cod2:&[char;3],ti:&mut[f64;3],tv:&mut[f64;3],l:&mut[f64;3],aa:&[usize;64],rl:&[Vec<f64>],pos_diff_flags:&[i32;3]) {
    let mut diff_indices = [0; 2]; let mut current_diff_count = 0;
    for k in 0..3 { if pos_diff_flags[k] == 1 { if current_diff_count < 2 { diff_indices[current_diff_count] = k; } current_diff_count += 1; } }
    if current_diff_count != 2 { return; }
    let d1_idx = diff_indices[0]; let d2_idx = diff_indices[1];
    let mut codint1 = *cod1; codint1[d1_idx] = cod2[d1_idx];
    let mut codint2 = *cod1; codint2[d2_idx] = cod2[d2_idx];
    let aa1=aa[num(cod1)]; let aa2=aa[num(cod2)];
    let aaint1_path1=aa[num(&codint1)]; let aaint1_path2=aa[num(&codint2)];
    let l_path1 = rl[aa1][aaint1_path1] * rl[aaint1_path1][aa2]; let l_path2 = rl[aa1][aaint1_path2] * rl[aaint1_path2][aa2];
    let p1 = if (l_path1 + l_path2) != 0.0 { l_path1 / (l_path1 + l_path2) } else { 0.5 }; let p2 = 1.0 - p1;
    let mut non_diff_pos_plus_1 = 0;
    for k_idx in 0..3 { if pos_diff_flags[k_idx] == 0 { non_diff_pos_plus_1 = k_idx + 1; break; } }
    if non_diff_pos_plus_1 > 0 {
        l[catsite(cod1[0],cod1[1],cod1[2],non_diff_pos_plus_1 as i32)]+=0.333333; l[catsite(cod2[0],cod2[1],cod2[2],non_diff_pos_plus_1 as i32)]+=0.333333;
        l[catsite(codint1[0],codint1[1],codint1[2],non_diff_pos_plus_1 as i32)]+=0.333333*p1; l[catsite(codint2[0],codint2[1],codint2[2],non_diff_pos_plus_1 as i32)]+=0.333333*p2;
    }
    titv1(cod1,&codint1,p1,ti,tv,l); titv1(&codint1,cod2,p1,ti,tv,l);
    titv1(cod1,&codint2,p2,ti,tv,l); titv1(&codint2,cod2,p2,ti,tv,l);
}
fn titv3(cod1:&[char;3],cod2:&[char;3],ti:&mut[f64;3],tv:&mut[f64;3],l:&mut[f64;3],aa:&[usize;64],rl:&[Vec<f64>]) {
    let mut codint1_paths = [['A'; 3]; 6]; let mut codint2_paths = [['A'; 3]; 6];
    let mut like = [0.0; 6]; let mut path_idx = 0;
    for i in 0..3 { for j in 0..3 {
        if j == i { continue; }
        let mut c2_intermediate = *cod1; let mut c3_intermediate;
        c2_intermediate[i] = cod2[i]; codint1_paths[path_idx] = c2_intermediate;
        c3_intermediate = c2_intermediate; c3_intermediate[j] = cod2[j]; codint2_paths[path_idx] = c3_intermediate;
        let aa_orig = aa[num(cod1)]; let aa_c2_int = aa[num(&c2_intermediate)];
        let aa_c3_int = aa[num(&c3_intermediate)]; let aa_final = aa[num(cod2)];
        like[path_idx] = rl[aa_orig][aa_c2_int] * rl[aa_c2_int][aa_c3_int] * rl[aa_c3_int][aa_final]; path_idx += 1;
    } }
    let somli:f64 = like.iter().sum(); let mut p = [0.0; 6];
    if somli > 0.0 { for i in 0..6 { p[i] = like[i] / somli; } } else { for i in 0..6 { p[i] = 1.0 / 6.0; } }
    for i in 0..6 { if p[i] > 0.0 { titv1(cod1, &codint1_paths[i], p[i], ti, tv, l); titv1(&codint1_paths[i], &codint2_paths[i], p[i], ti, tv, l); titv1(&codint2_paths[i], cod2, p[i], ti, tv, l); } }
}
fn prefastlwl(rl:&[Vec<f64>], tl0:&mut Vec<Vec<f64>>, tl1:&mut Vec<Vec<f64>>, tl2:&mut Vec<Vec<f64>>, tti0:&mut Vec<Vec<f64>>, tti1:&mut Vec<Vec<f64>>, tti2:&mut Vec<Vec<f64>>, ttv0:&mut Vec<Vec<f64>>, ttv1:&mut Vec<Vec<f64>>, ttv2:&mut Vec<Vec<f64>>) {
    let mut aa_li_map=[0;64]; fill_aa(&mut aa_li_map);
    for i in 0..64 { for j in i..64 {
        let mut l_accum = [0.0;3]; let mut ti_accum = [0.0;3]; let mut tv_accum = [0.0;3];
        let cod1_chars = decode_codon(i); let cod2_chars = decode_codon(j);
        let mut pos_diff_flags=[0i32;3]; let mut nbdiff=0;
        for p in 0..3 { if cod1_chars[p]!=cod2_chars[p] { nbdiff+=1; pos_diff_flags[p]=1; } }
        if nbdiff == 0 { for p_idx in 0..3 { let site_type = catsite(cod1_chars[0], cod1_chars[1], cod1_chars[2], (p_idx + 1) as i32); l_accum[site_type] += 1.0; }
        } else if nbdiff == 1 {
            for p_idx in 0..3 { if pos_diff_flags[p_idx] == 0 { let site_type1 = catsite(cod1_chars[0], cod1_chars[1], cod1_chars[2], (p_idx + 1) as i32); l_accum[site_type1] += 0.5; let site_type2 = catsite(cod2_chars[0], cod2_chars[1], cod2_chars[2], (p_idx + 1) as i32); l_accum[site_type2] += 0.5; } }
            titv1(&cod1_chars,&cod2_chars,1.0,&mut ti_accum,&mut tv_accum,&mut l_accum);
        } else if nbdiff==2 {
            for p_idx in 0..3 { if pos_diff_flags[p_idx] == 0 { let site_type1 = catsite(cod1_chars[0], cod1_chars[1], cod1_chars[2], (p_idx + 1) as i32); l_accum[site_type1] += 0.5; let site_type2 = catsite(cod2_chars[0], cod2_chars[1], cod2_chars[2], (p_idx + 1) as i32); l_accum[site_type2] += 0.5; break; } }
            titv2(&cod1_chars,&cod2_chars,&mut ti_accum,&mut tv_accum,&mut l_accum,&aa_li_map,rl,&pos_diff_flags);
        } else if nbdiff==3 { titv3(&cod1_chars,&cod2_chars,&mut ti_accum,&mut tv_accum,&mut l_accum,&aa_li_map,rl); }
        tl0[i][j]=l_accum[0]; tl0[j][i]=l_accum[0]; tl1[i][j]=l_accum[1]; tl1[j][i]=l_accum[1]; tl2[i][j]=l_accum[2]; tl2[j][i]=l_accum[2];
        tti0[i][j]=ti_accum[0]; tti0[j][i]=ti_accum[0]; tti1[i][j]=ti_accum[1]; tti1[j][i]=ti_accum[1]; tti2[i][j]=ti_accum[2]; tti2[j][i]=ti_accum[2];
        ttv0[i][j]=tv_accum[0]; ttv0[j][i]=tv_accum[0]; ttv1[i][j]=tv_accum[1]; ttv1[j][i]=tv_accum[1]; ttv2[i][j]=tv_accum[2]; ttv2[j][i]=tv_accum[2];
    } }
}
fn fastlwl_pair(codon_indices1:&[u8], codon_indices2:&[u8], tti0:&[Vec<f64>], tti1:&[Vec<f64>], tti2:&[Vec<f64>], ttv0:&[Vec<f64>], ttv1:&[Vec<f64>], ttv2:&[Vec<f64>], tl0:&[Vec<f64>], tl1:&[Vec<f64>], tl2:&[Vec<f64>]) -> (f64,f64) {
    let nb_codon = codon_indices1.len(); let mut l_sum=[0.0;3]; let mut ti_sum=[0.0;3]; let mut tv_sum=[0.0;3];
    for k in 0..nb_codon {
        let num1 = codon_indices1[k] as usize; let num2 = codon_indices2[k] as usize;
        if num1==64||num2==64 {continue;}
        l_sum[0]+=tl0[num1][num2]; l_sum[1]+=tl1[num1][num2]; l_sum[2]+=tl2[num1][num2];
        ti_sum[0]+=tti0[num1][num2]; ti_sum[1]+=tti1[num1][num2]; ti_sum[2]+=tti2[num1][num2];
        tv_sum[0]+=ttv0[num1][num2]; tv_sum[1]+=ttv1[num1][num2]; tv_sum[2]+=ttv2[num1][num2];
    }
    let mut p_k=[0.0;3]; let mut q_k=[0.0;3];
    for k_type in 0..3 { if l_sum[k_type] == 0.0 { p_k[k_type] = 0.0; q_k[k_type] = 0.0; } else { p_k[k_type]=ti_sum[k_type]/l_sum[k_type]; q_k[k_type]=tv_sum[k_type]/l_sum[k_type]; } }
    let mut a_k_val=[0.0;3]; let mut b_k_val=[0.0;3];
    for k_type in 0..3 {
        if l_sum[k_type] == 0.0 { a_k_val[k_type] = f64::NAN; b_k_val[k_type] = f64::NAN; continue; }
        let denom_a = 1.0 - 2.0 * p_k[k_type] - q_k[k_type]; let denom_b = 1.0 - 2.0 * q_k[k_type];
        if denom_a <= 1e-15 || denom_b <= 1e-15 { 
            a_k_val[k_type] = f64::NAN; b_k_val[k_type] = f64::NAN;
        } else {
            a_k_val[k_type] = -0.5 * denom_a.ln() - 0.25 * denom_b.ln();
            b_k_val[k_type] = -0.25 * denom_b.ln();
        }
    }
    let mut ks_val = f64::NAN; let mut ka_val = f64::NAN;
    let l0_s = l_sum[0]; let l2_s = l_sum[1]; let l4_s = l_sum[2];
    let a0_li = a_k_val[0]; let a2_li = a_k_val[1]; let a4_li = a_k_val[2];
    let b0_li = b_k_val[0]; let b2_li = b_k_val[1]; let b4_li = b_k_val[2];
    
    // **FIXED (Sugerencia #2)**: Cálculo robusto de Ka
    if l0_s > 0.0 && !a0_li.is_nan() && !b0_li.is_nan() {
        if l2_s > 0.0 && !b2_li.is_nan() {
            ka_val = a0_li + (l0_s * b0_li + l2_s * b2_li) / (l0_s + l2_s);
        } else {
            ka_val = a0_li + b0_li;
        }
    }
    
    // **FIXED (Sugerencia #1)**: Cálculo robusto de Ks
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

    // Reajuste posterior para evitar valores negativos por ruido numérico
    if ka_val.is_finite() && ka_val < 0.0 { ka_val = 0.0; }
    if ks_val.is_finite() && ks_val < 0.0 { ks_val = 0.0; }

    (ka_val, ks_val)
}
fn get_li_values(codon_indices1:&[u8], codon_indices2:&[u8], tti0:&[Vec<f64>], tti1:&[Vec<f64>], tti2:&[Vec<f64>], ttv0:&[Vec<f64>], ttv1:&[Vec<f64>], ttv2:&[Vec<f64>], tl0:&[Vec<f64>], tl1:&[Vec<f64>], tl2:&[Vec<f64>]) -> (f64, f64) {
    fastlwl_pair(codon_indices1, codon_indices2, tti0, tti1, tti2, ttv0, ttv1, ttv2, tl0, tl1, tl2)
}
fn calculate_syn_nonsyn_from_indices(codon_indices1: &[u8], codon_indices2: &[u8]) -> (f64, f64) {
    let mut count_valid_codons = 0; let mut syn_diffs = 0.0; let mut nonsyn_diffs = 0.0;
    let mut sum_syn_sites_seq1 = 0.0; let mut sum_syn_sites_seq2 = 0.0;
    for k in 0..codon_indices1.len() {
        let idx1 = codon_indices1[k] as usize; let idx2 = codon_indices2[k] as usize;
        if idx1 >= 64 || idx2 >= 64 { continue; }
        count_valid_codons += 1; sum_syn_sites_seq1 += SYN_SITE_ARRAY[idx1] as f64; sum_syn_sites_seq2 += SYN_SITE_ARRAY[idx2] as f64;
        if idx1 != idx2 { let aa1 = AA_ARRAY[idx1]; let aa2 = AA_ARRAY[idx2]; if aa1 == aa2 { syn_diffs += 1.0; } else { nonsyn_diffs += 1.0; } }
    }
    if count_valid_codons == 0 { return (f64::NAN, f64::NAN); }
    let potential_syn_sites = (sum_syn_sites_seq1/3.0 + sum_syn_sites_seq2/3.0) / 2.0;
    let potential_nonsyn_sites = (count_valid_codons as f64) * 3.0 - potential_syn_sites;
    let ps = if potential_syn_sites > 0.0 { syn_diffs / potential_syn_sites } else { 0.0 };
    let pn = if potential_nonsyn_sites > 0.0 { nonsyn_diffs / potential_nonsyn_sites } else { 0.0 };
    
    let mut ds = if ps < 0.749 { -0.75 * (1.0 - (4.0/3.0) * ps).ln() } else { f64::NAN };
    let mut dn = if pn < 0.749 { -0.75 * (1.0 - (4.0/3.0) * pn).ln() } else { f64::NAN };

    // Reajuste posterior para evitar valores negativos
    if ds.is_finite() && ds < 0.0 { ds = 0.0; }
    if dn.is_finite() && dn < 0.0 { dn = 0.0; }

    (dn, ds)
}

// --- NUEVO: Funciones para la codificación compacta de ADN ---

/// Convierte un byte de base (ASCII) a un número compacto (0-4).
/// A=0, C=1, G=2, T/U=3, Otro=4 (Ambiguo/N)
#[inline(always)]
fn base_to_dna5(b: u8) -> u8 {
    const LUT: [u8; 256] = {
        let mut t = [4; 256]; // 4 es el valor para N o ambiguos
        t[b'A' as usize] = 0; t[b'a' as usize] = 0;
        t[b'C' as usize] = 1; t[b'c' as usize] = 1;
        t[b'G' as usize] = 2; t[b'g' as usize] = 2;
        t[b'T' as usize] = 3; t[b't' as usize] = 3;
        t[b'U' as usize] = 3; t[b'u' as usize] = 3;
        t
    };
    LUT[b as usize]
}

/// Convierte una secuencia de bytes ASCII a una secuencia compacta Dna5.
fn seq_to_dna5(bytes: &[u8]) -> Vec<u8> {
    bytes.iter().map(|&b| base_to_dna5(b)).collect()
}

/// Convierte una secuencia compacta Dna5 a una lista de índices de codones.
fn dna5_to_codon_indices(dna5_bases: &[u8], model: Model) -> Vec<u8> {
    let mut v = Vec::with_capacity(dna5_bases.len() / 3);
    for chunk in dna5_bases.chunks_exact(3) {
        let (b1, b2, b3) = (chunk[0], chunk[1], chunk[2]);

        if b1 > 3 || b2 > 3 || b3 > 3 { // Si alguna base es ambigua
            v.push(64);
        } else {
            let index = match model {
                // Modelo Li usa orden A,C,G,T (0,1,2,3), que coincide con Dna5
                Model::Li => 16 * b1 + 4 * b2 + b3,
                // **FIXED**: Modelo Nei usa T(0),C(1),A(2),G(3).
                // Se remapea desde nuestro Dna5 (A=0,C=1,G=2,T=3) al orden de Nei.
                // Y se usa la fórmula de indexación estándar.
                Model::Nei => {
                    let remap = |b: u8| match b { 0 => 2, 1 => 1, 2 => 3, 3 => 0, _ => unreachable!() };
                    let n1 = remap(b1); // Base en pos 0
                    let n2 = remap(b2); // Base en pos 1
                    let n3 = remap(b3); // Base en pos 2
                    16 * n1 + 4 * n2 + n3
                }
            };
            v.push(index);
        }
    }
    v
}


#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
enum Model {
    Nei,
    Li,
}

#[derive(Clone, Copy, Debug)]
struct DsDn {
    dn: f64,
    ds: f64,
}

/// Calcula dN/dS para secuencias usando los modelos Nei-Gojobori o Li (1993) de forma eficiente.
#[derive(Parser, Debug)]
#[command(version = "1.0.0", about, long_about = None)]
struct Args {
    /// Archivo de entrada con secuencias en formato FASTA
    #[arg(required = true)]
    input_file: String,

    /// Nombre base para los archivos de salida
    #[arg(short, long, default_value = "output")]
    output: String,

    /// Número de procesos paralelos
    #[arg(short, long, default_value_t = 4)]
    workers: usize,

    /// Si se especifica, calcula la media de dN y dS agrupada por linaje contra todos los demás
    #[arg(long, group = "output_mode")]
    lineage: bool,

    /// Si se especifica, calcula la media de dN/dS entre grupos predefinidos
    #[arg(long, group = "output_mode")]
    group_average: bool,

    /// En el modo linaje, agrupa por la primera letra del ID de la secuencia (ignora la división por '_')
    #[arg(long, requires = "lineage")]
    first_letter_lineage: bool,

    /// Modelo a utilizar
    #[arg(long, value_enum, default_value_t = Model::Nei)]
    model: Model,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let args = Args::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.workers)
        .stack_size(4 * 1024 * 1024)
        .build_global()?;

    // --- Lectura y Deduplicación con codificación Dna5 ---
    info!("Leyendo y deduplicando secuencias desde: {}", args.input_file);
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut seq_to_ids: FxHashMap<Vec<u8>, Vec<Vec<u8>>> = FxHashMap::default();
    let mut reader = parse_fastx_file(&args.input_file)?;
    while let Some(record) = reader.next() {
        let rec = record?;
        let seq_dna5 = seq_to_dna5(&rec.seq()); // Convertir a Dna5 al momento
        let id = rec.id().to_vec();
        ids.push(id.clone());
        seq_to_ids.entry(seq_dna5).or_default().push(id);
    }
    info!("Se encontraron {} secuencias en total.", ids.len());
    if ids.is_empty() {
        eprintln!("No se encontraron secuencias en el archivo de entrada.");
        std::process::exit(1);
    }
    
    // Ahora las claves del mapa (secuencias únicas) están en formato Dna5, ahorrando memoria.
    let mut unique_seqs_dna5: Vec<Vec<u8>> = seq_to_ids.keys().cloned().collect();
    unique_seqs_dna5.sort_unstable(); 
    let n_u = unique_seqs_dna5.len();
    info!("Se encontraron {} secuencias únicas.", n_u);

    // --- Mapeo de IDs a índices únicos ---
    let id_to_uidx: FxHashMap<Vec<u8>, usize> = unique_seqs_dna5.iter().enumerate()
        .flat_map(|(u_idx, seq_dna5)| {
            seq_to_ids.get(seq_dna5).unwrap().iter().map(move |id| (id.clone(), u_idx))
        })
        .collect();
    
    // --- Codificación de Secuencias Únicas ---
    info!("Codificando secuencias únicas a índices de codones...");
    let unique_codon_indices: Arc<Vec<Vec<u8>>> = Arc::new(unique_seqs_dna5
        .par_iter()
        .map(|s_dna5| dna5_to_codon_indices(s_dna5, args.model))
        .collect());
    drop(unique_seqs_dna5); 
    drop(seq_to_ids); 

    // --- Precomputación para Modelo Li ---
    let (tl0, tl1, tl2, tti0, tti1, tti2, ttv0, ttv1, ttv2);
    if args.model == Model::Li {
        info!("Precalculando tablas para el modelo Li (1993)...");
        let mut rl_li = vec![vec![0.0; 21]; 21];
        for i in 2..=20 { for j in 1..i { rl_li[i][j] = MAT[j-1][i-2]; } }
        for i in 1..=20 { rl_li[i][i] = 1.0; for j in i+1..=20 { rl_li[i][j] = rl_li[j][i]; } }
        let mut minrl = rl_li.get(1).and_then(|row| row.get(1)).copied().unwrap_or(0.01);
        for i in 1..=20 { for j in i+1..=20 { if rl_li[i][j]<minrl {minrl=rl_li[i][j];} } }
        for i in 0..=20 { rl_li[0][i]=minrl; rl_li[i][0]=minrl; }
        
        let mut tl0_mut = vec![vec![0.0;64];64]; let mut tl1_mut = vec![vec![0.0;64];64]; let mut tl2_mut = vec![vec![0.0;64];64];
        let mut tti0_mut = vec![vec![0.0;64];64]; let mut tti1_mut = vec![vec![0.0;64];64]; let mut tti2_mut = vec![vec![0.0;64];64];
        let mut ttv0_mut = vec![vec![0.0;64];64]; let mut ttv1_mut = vec![vec![0.0;64];64]; let mut ttv2_mut = vec![vec![0.0;64];64];
        prefastlwl(&rl_li, &mut tl0_mut, &mut tl1_mut, &mut tl2_mut, &mut tti0_mut, &mut tti1_mut, &mut tti2_mut, &mut ttv0_mut, &mut ttv1_mut, &mut ttv2_mut);
        
        tl0 = Arc::new(tl0_mut); tl1 = Arc::new(tl1_mut); tl2 = Arc::new(tl2_mut); 
        tti0 = Arc::new(tti0_mut); tti1 = Arc::new(tti1_mut); tti2 = Arc::new(tti2_mut); 
        ttv0 = Arc::new(ttv0_mut); ttv1 = Arc::new(ttv1_mut); ttv2 = Arc::new(ttv2_mut);
        info!("Precalculación para Li finalizada.");
    } else {
        tl0 = Arc::new(Vec::new()); tl1 = Arc::new(Vec::new()); tl2 = Arc::new(Vec::new()); 
        tti0 = Arc::new(Vec::new()); tti1 = Arc::new(Vec::new()); tti2 = Arc::new(Vec::new()); 
        ttv0 = Arc::new(Vec::new()); ttv1 = Arc::new(Vec::new()); ttv2 = Arc::new(Vec::new());
    }

    // --- FASE 1: Cálculo en paralelo de pares únicos en memoria ---
    let total_unique_pairs = if n_u > 0 { n_u * (n_u - 1) / 2 } else { 0 };
    info!("Calculando dN/dS para {} pares de secuencias únicas...", total_unique_pairs);
    let mut pair_results: Vec<DsDn> = vec![DsDn { dn: f64::NAN, ds: f64::NAN }; total_unique_pairs];
    
    let pb_calc = ProgressBar::new(n_u as u64);
    pb_calc.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Calculando pares únicos: {pos}/{len} ({eta})")
        .progress_chars("#>-"));

    let pair_results_ptr_addr = pair_results.as_mut_ptr() as usize;
    (0..n_u).into_par_iter().for_each(|i| {
        let pair_results_ptr = pair_results_ptr_addr as *mut DsDn;
        let cod_i = &unique_codon_indices[i];

        for j in (i + 1)..n_u {
            let cod_j = &unique_codon_indices[j];
            let (dn, ds) = if args.model == Model::Li {
                get_li_values(cod_i, cod_j, &tti0, &tti1, &tti2, &ttv0, &ttv1, &ttv2, &tl0, &tl1, &tl2)
            } else {
                calculate_syn_nonsyn_from_indices(cod_i, cod_j)
            };
            
            let flat_idx = (i * n_u + j) - ((i + 1) * (i + 2)) / 2;
            unsafe {
                *pair_results_ptr.add(flat_idx) = DsDn { dn, ds };
            }
        }
        pb_calc.inc(1);
    });
    pb_calc.finish_with_message("Cálculo de pares únicos completado.");
    
    // --- FASE 2: Generación de la salida a partir de los resultados en memoria ---
    let get_result = |u_i: usize, u_j: usize| -> DsDn {
        if u_i == u_j {
            DsDn { dn: 0.0, ds: 0.0 }
        } else {
            let (i, j) = if u_i < u_j { (u_i, u_j) } else { (u_j, u_i) };
            let flat_idx = (i * n_u + j) - ((i + 1) * (i + 2)) / 2;
            pair_results[flat_idx]
        }
    };
    
    if args.group_average {
        info!("Calculando promedios de grupo dN/dS...");
        let mut group_map: FxHashMap<Vec<u8>, Vec<Vec<u8>>> = FxHashMap::default();
        for id_bytes in &ids {
            let group_key = id_bytes.split(|&b| b == b'_').next().unwrap_or(id_bytes).to_vec();
            group_map.entry(group_key).or_default().push(id_bytes.clone());
        }

        let group_keys: Vec<Vec<u8>> = group_map.keys().cloned().collect();
        let mut group_pairs_to_process = Vec::new();
        for (i, g1_key) in group_keys.iter().enumerate() {
            for g2_key in group_keys.iter().skip(i) {
                group_pairs_to_process.push((g1_key.clone(), g2_key.clone()));
            }
        }

        let pb_group = ProgressBar::new(group_pairs_to_process.len() as u64);
        pb_group.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Grupos: {pos}/{len} ({eta})")
                .progress_chars("#>-"),
        );

        let group_results: Vec<String> = group_pairs_to_process
            .par_iter()
            .map(|(g1_key, g2_key)| {
                let ids1 = &group_map[g1_key];
                let ids2 = &group_map[g2_key];
                let mut pair_dn_ds_ratios = Vec::new();

                for (i, id1) in ids1.iter().enumerate() {
                    let u_idx1 = id_to_uidx[id1];
                    
                    let iter_ids2: Box<dyn Iterator<Item = &Vec<u8>>> = if g1_key == g2_key {
                        Box::new(ids2.iter().skip(i + 1))
                    } else {
                        Box::new(ids2.iter())
                    };

                    for id2 in iter_ids2 {
                        let u_idx2 = id_to_uidx[id2];
                        let result = get_result(u_idx1, u_idx2);
                        
                        if result.dn.is_finite() && result.ds.is_finite() && result.ds > 0.0 {
                            pair_dn_ds_ratios.push(result.dn / result.ds);
                        }
                    }
                }
                
                pb_group.inc(1);
                if pair_dn_ds_ratios.is_empty() {
                    format!("{}\t{}\t{}\t{}\t0\tNaN\tNaN\t[NaN, NaN]\n", String::from_utf8_lossy(g1_key), String::from_utf8_lossy(g2_key), ids1.len(), ids2.len())
                } else {
                    let num_comparisons = pair_dn_ds_ratios.len();
                    let mean_dn_ds: f64 = pair_dn_ds_ratios.iter().sum::<f64>() / num_comparisons as f64;
                    let variance: f64 = if num_comparisons > 1 {
                        pair_dn_ds_ratios.iter().map(|&val| (val - mean_dn_ds).powi(2)).sum::<f64>() / (num_comparisons - 1) as f64
                    } else { 0.0 };
                    let standard_error = if num_comparisons > 0 {(variance / num_comparisons as f64).sqrt()} else {0.0};
                    let t_value = 1.96; 
                    let ci_half_width = t_value * standard_error;
                    let ci_lower = mean_dn_ds - ci_half_width;
                    let ci_upper = mean_dn_ds + ci_half_width;
                    format!("{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t[{:.6}, {:.6}]\n",
                        String::from_utf8_lossy(g1_key), String::from_utf8_lossy(g2_key), ids1.len(), ids2.len(), num_comparisons,
                        mean_dn_ds, standard_error, ci_lower, ci_upper
                    )
                }
            }).collect();

        let output_path_group = format!("{}_group_avg_dn_ds.tsv", args.output);
        let mut out_file_group = BufWriter::new(File::create(&output_path_group)?);
        writeln!(out_file_group, "Group1\tGroup2\tNumSeqs1\tNumSeqs2\tNumComparisons\tMean_dN/dS\tStdError\t95%CI")?;
        for line in group_results {
            if !line.is_empty() { out_file_group.write_all(line.as_bytes())?; }
        }
        pb_group.finish_with_message("Cálculo de promedios de grupo completado.");
        info!("Resultados guardados en {}", output_path_group);

    } else if args.lineage {
        info!("Calculando resumen de dN/dS por linaje...");
        
        let lineage_results: Vec<String> = (0..ids.len())
            .into_par_iter()
            .map(|i| {
                let id_i_bytes = &ids[i];
                let u_i = id_to_uidx[id_i_bytes];
                
                // Mapa para este genoma: Linaje de comparación -> (Suma dn, Suma ds, Conteo)
                let mut local_aggr: FxHashMap<Vec<u8>, (f64, f64, usize)> = FxHashMap::default();

                for j in 0..ids.len() {
                    if i == j { continue; }
                    let id_j_bytes = &ids[j];
                    let u_j = id_to_uidx[id_j_bytes];
                    
                    let result = get_result(u_i, u_j);

                    if result.dn.is_finite() && result.ds.is_finite() {
                        let lineage_key = if args.first_letter_lineage {
                            id_j_bytes.iter().next().map(|&b| vec![b]).unwrap_or_default()
                        } else {
                            id_j_bytes.split(|&b| b == b'_').next().unwrap_or(id_j_bytes).to_vec()
                        };
                        
                        let entry = local_aggr.entry(lineage_key).or_default();
                        entry.0 += result.dn;
                        entry.1 += result.ds;
                        entry.2 += 1;
                    }
                }

                let mut output_lines = String::new();
                for (lineage_key, (sum_dn, sum_ds, count)) in local_aggr {
                    let mean_dn = sum_dn / count as f64;
                    let mean_ds = sum_ds / count as f64;
                    let ratio = if mean_ds == 0.0 {
                        if mean_dn == 0.0 { 0.0 } else { f64::INFINITY }
                    } else {
                        mean_dn / mean_ds
                    };
                    output_lines.push_str(&format!("{}\t{}\t{:.6}\t{:.6}\t{:.6}\n", 
                        String::from_utf8_lossy(id_i_bytes),
                        String::from_utf8_lossy(&lineage_key),
                        mean_dn, mean_ds, ratio));
                }
                output_lines
            })
            .collect();
            
        let output_path_lineage = format!("{}_lineage_summary.tsv", args.output);
        let mut out_file_lineage = BufWriter::new(File::create(&output_path_lineage)?);
        writeln!(out_file_lineage, "Genome\tAgainst_Lineage\tMean_dN\tMean_dS\tdN/dS_Ratio")?;
        for block in lineage_results {
            out_file_lineage.write_all(block.as_bytes())?;
        }
        info!("Resumen de linaje guardado en {}", output_path_lineage);

    } else {
        // Modo por defecto: pairwise
        info!("Generando informe de resultados por pares...");
        let total_pairs_to_write = ids.len() * (ids.len() - 1) / 2;
        let pb_write = ProgressBar::new(total_pairs_to_write as u64);
        pb_write.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] Escribiendo pares: {pos}/{len} ({eta})")
            .progress_chars("#>-"));

        let (tx, rx) = unbounded::<String>();

        let writer_thread = thread::spawn({
            let output_path = format!("{}_pairwise_results.tsv", args.output);
            let model = args.model;
            move || -> Result<(), std::io::Error> {
                let mut out_file = BufWriter::new(File::create(output_path)?);
                let header = match model {
                    Model::Li => "Seq1\tSeq2\tdN(Ka)\tdS(Ks)\tdN/dS\n",
                    Model::Nei => "Seq1\tSeq2\tdN\tdS\tdN/dS\n",
                };
                out_file.write_all(header.as_bytes())?;
                for line_batch in rx {
                    out_file.write_all(line_batch.as_bytes())?;
                }
                Ok(())
            }
        });

        (0..ids.len()).into_par_iter().for_each_with(tx, |s, i| {
            let mut local_buffer = String::with_capacity(1024 * 8);
            let id_i_bytes = &ids[i];
            let u_i = id_to_uidx[id_i_bytes];

            for j in (i + 1)..ids.len() {
                let id_j_bytes = &ids[j];
                let u_j = id_to_uidx[id_j_bytes];
                
                let result = get_result(u_i, u_j);
                let ratio = if result.ds == 0.0 {
                    if result.dn == 0.0 { 0.0 } else { f64::INFINITY }
                } else {
                    result.dn / result.ds
                };
                
                let line = format!("{}\t{}\t{:.6}\t{:.6}\t{:.6}\n",
                    String::from_utf8_lossy(id_i_bytes),
                    String::from_utf8_lossy(id_j_bytes),
                    result.dn, result.ds, ratio);
                local_buffer.push_str(&line);

                if local_buffer.len() > 1024 * 4 {
                    s.send(local_buffer.clone()).unwrap();
                    local_buffer.clear();
                }
            }
            if !local_buffer.is_empty() {
                s.send(local_buffer).unwrap();
            }
            pb_write.inc((ids.len() - 1 - i) as u64);
        });
        
        // `for_each_with` automatically drops the sender `tx` when it's done.
        
        writer_thread.join().unwrap().unwrap();
        pb_write.finish_with_message("Escritura de pares completada.");
        info!("Resultados guardados en {}_pairwise_results.tsv", args.output);
    }

    info!("Todos los procesos han finalizado con éxito.");
    Ok(())
}
