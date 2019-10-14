import os
from nest.parsers.vcfReader import Reader 



class Writer:

    def __init__(self, out_path):
        self.out_path = out_path
        self.vcf_writer = open(out_path, 'w')

    def writeHeaders(self, headers):
        for header in headers:
            if header == 'fileFormat':
                self.vcf_writer.write('##fileformat={0}\n'.format(
                                                headers[header]))
            elif header == 'info':
                for info in headers['info']:
                    if info.source == None:
                        self.vcf_writer.write(('##INFO=<ID={0},'.format(
                                                                    info.id)
                                + 'Number={0},'.format(info.number)
                                + 'Type={0},'.format(info.type)
                                + 'Description="{0}">\n'.format(
                                                        info.description)
                                ))
                    else:
                        self.vcf_writer.write(('##INFO=<ID={0},'.format(
                                                                    info.id)
                                + 'Number={0},'.format(info.number)
                                + 'Type={0},'.format(info.type)
                                + 'Description="{0}",'.format(info.description)
                                + 'Source="{0}",'.format(info.source)
                                + 'Version="{0}">\n'.format(info.version)))
            elif header == 'filter':
                for filters in headers['filter']:
                    self.vcf_writer.write(('##FILTER=<ID={0},'.format(
                                                            filters.id)
                        + 'Description="{0}">\n'.format(filters.description)
                        ))
            elif header == 'format':
                    for formats in headers['format']:
                        #print(formats.id)
                        self.vcf_writer.write(('##FORMAT=<ID={0},'.format(
                                                            formats.id)
                        + 'Number={0},'.format(formats.number)
                        + 'Type={0},'.format(formats.type)
                        + 'Description="{0}">\n'.format(formats.description)
                        ))
            elif header == 'alt':
                for alt in headers['alt']:
                    self.vcf_writer.write(('##ALT=<ID={0},'.format(alt.id)
                        + 'Description="{0}">\n'.format(alt.description)))
            elif header == 'assembly':
                self.vcf_writer.write('##assembly={0}\n'.format(
                                    headers[header]))
            elif header == 'contig':
                for contig in headers['contig']:
                    if contig.url == None:
                        self.vcf_writer.write(('##contig=<ID={0},'.format(
                                                            contig.id)
                                + 'length={0}>\n'.format(contig.length)
                                ))
                    else:
                        self.vcf_writer.write(('##contig=<ID={0},'.format(
                                                            contig.id)
                                + 'length={0},'.format(contig.length)
                                + 'URL={0}>\n'.format(contig.url)))
            elif header == 'random':
                for rands in headers['random']:
                    self.vcf_writer.write(('##{0}={1}\n'.format(
                                rands.field, rands.value)))
                
        ##Write header line
        header_line = headers['fields'] + headers['samples']
        self.vcf_writer.write('#{0}\n'.format('\t'.join(header_line)))

    def writeRecords(self, records):
        rec = list()
        #print(records.Samples)
        rec.append('{0}'.format(records.CHROM))
        rec.append('{0}'.format(records.POS))
        rec.append('{0}'.format(records.UID))
        rec.append(','.join(records.REF))
        alt = ','.join(records.ALT)
        rec.append(alt)
        rec.append('{0}'.format(records.QUAL))
        rec.append(records.FILTER)
        info = list()
        for key, value in records.INFO.items():
            #print(key,value)
            # If values are absent, do not write value
            if len(value) == 1:
                if value[0] is None:
                    continue
                elif type(value[0]) == bool and value[0]:
                    info.append('{0}'.format(key))
                elif type(value[0]) == bool and value[0]:
                    continue
                else:
                    info.append('{0}={1}'.format(key, value[0]))
            else:
                info.append('{0}={1}'.format(key, 
                            ','.join([str(val) for val in value])))
        rec.append(';'.join(info))
        all_formats = list()
        #print(records.CHROM, records.POS)
        #print(records.Samples)
        for key, values in records.Samples.items():
            all_formats = list(values.keys())
            break
        #rec.append(':'.join(formats))
        formats = list()
        for sample, info in records.Samples.items():
            #1/1:0,41:41:99:1535,123,0
            #1/1:1:0,41:0:41:4:99:9:1535,123,
            sam_string = list()
            for keys in all_formats:
                value = info[keys]
                if len(value) == 1:
                    sam_string.append(str(value[0]))
                else:
                    sstring = [str(val) for val in value]
                    sam_string.append(','.join(sstring))
                if key not in formats:
                    formats.append(keys)
                #print(value)
                #if value[0] == None:
                #    continue
                #elif type(value[0]) == list and None in value[0]:
                #    continue
                #elif type(value[0]) == list and len(value[0]) == 1:
                #    sam_string.append(str(value[0][0]))
                #    if keys not in formats:
                #        formats.append(keys)
                #elif type(value[0]) == list and len(value[0]) > 1 :
                #    value = ','.join([
                #    str(val) if val != None else 'NaN'for val in value[0]])
                #    sam_string.append(value)
                #    if keys not in formats:
                #        formats.append(keys)
                #else:
                #    if value[0] == None:
                #        continue
                #    else:
                #        sam_string.append(str(value[0]))
                #        if keys not in formats:
                #            formats.append(keys)
            #print(':'.join(formats))
            #print(':'.join(sam_string))
            rec.append(':'.join(formats))
            rec.append(':'.join(sam_string))
        self.vcf_writer.write('{0}\n'.format('\t'.join(rec)))
        return

    def closeWriter(self):
        self.vcf_writer.close()
        return
