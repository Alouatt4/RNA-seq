import sys

#convert gff3 to gtf
def convert_edta_to_gtf(gff_in, gtf_out):
  with open(gff_in, 'r') as infile, open(gtf_out, 'w') as outfile:
    count = 0
    for line in infile:
      if line.startswith('#') or not line.strip():
        continue

      parts = line.strip().split('\t')
      if len(parts) < 9:
        continue

      chrom = parts[0]
      source = 'EDTA'
      feature = 'exon'  # TEtranscripts exige "exon" na 3ª coluna
      start = parts[3]
      end = parts[4]
      score = parts[5] if parts[5] != '.' else '0'
      strand = parts[6] if parts[6] in ['+', '-'] else '.'
      frame = '.'

      # Parse da coluna 9 do GFF3
      attr_list = parts[8].split(';')
      attrs = {}
      for item in attr_list:
        if '=' in item:
          k, v = item.split('=', 1)
          attrs[k] = v

      te_name = attrs.get('Name', f'TE_{count}')
      classification = attrs.get('Classification', 'Unknown')

      # Separa Classe e Família (ex: LTR/Copia -> class_id=LTR, family_id=Copia)
      if '/' in classification:
        te_class, te_family = classification.split('/', 1)
      else:
        te_class = classification
        te_family = classification

      count += 1
      transcript_id = f'{te_name}_dup{count}'
      gene_id = te_name

      gtf_attributes = (
          f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; family_id'
          f' "{te_family}"; class_id "{te_class}";'
      )

      outfile.write(
          f'{chrom}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{gtf_attributes}\n'
      )


if __name__ == '__main__':
  if len(sys.argv) != 3:
    print('Uso: python convert_edta_to_tetranscripts.py <input.gff3> <output.gtf>')
    sys.exit(1)
  convert_edta_to_gtf(sys.argv[1], sys.argv[2])
