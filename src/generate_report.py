#!/usr/bin/env python3
"""
Generate a project report for small RNA sequencing data
with metadata extracted from FASTQ files.
"""

import os
import glob
import gzip
import re
from datetime import datetime
from pathlib import Path
from collections import defaultdict
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak, Image
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT


class FASTQMetadataExtractor:
    """Extract metadata from FASTQ files"""
    
    def __init__(self, fastq_dir):
        self.fastq_dir = fastq_dir
        self.samples = {}
        self.parse_samples()
    
    def parse_samples(self):
        """Parse sample names and file sizes from FASTQ directory"""
        r1_files = sorted(glob.glob(os.path.join(self.fastq_dir, "*-R1.fastq.gz")))
        r2_files = sorted(glob.glob(os.path.join(self.fastq_dir, "*-R2.fastq.gz")))
        
        for r1_file in r1_files:
            basename = os.path.basename(r1_file).replace("-R1.fastq.gz", "")
            r2_file = os.path.join(self.fastq_dir, basename + "-R2.fastq.gz")
            
            r1_size = os.path.getsize(r1_file) / (1024**3)  # Convert to GB
            r2_size = os.path.getsize(r2_file) / (1024**3) if os.path.exists(r2_file) else 0
            
            # Extract i7 and i5 indices from basename
            # Format: xRXXX-LX-GX-PXX-XXXXXXXX-XXXXXXXX
            parts = basename.split('-')
            i7_index = parts[-2] if len(parts) >= 2 else 'N/A'
            i5_index = parts[-1] if len(parts) >= 1 else 'N/A'
            
            # Count reads in R1 file (sample first few lines)
            read_count = self._count_reads_sample(r1_file)
            
            self.samples[basename] = {
                'r1_path': r1_file,
                'r2_path': r2_file if os.path.exists(r2_file) else None,
                'r1_size_gb': r1_size,
                'r2_size_gb': r2_size,
                'read_count': read_count,
                'total_size_gb': r1_size + r2_size,
                'i7_index': i7_index,
                'i5_index': i5_index
            }
    
    def _count_reads_sample(self, fastq_file, sample_lines=1000):
        """Estimate read count by sampling lines from gzip FASTQ file"""
        try:
            with gzip.open(fastq_file, 'rt') as f:
                count = 0
                for _ in range(sample_lines):
                    line = f.readline()
                    if not line:
                        break
                    if line.startswith('@') and not line.startswith('@+'):
                        count += 1
            
            # Estimate total from sample
            total_size = os.path.getsize(fastq_file)
            file_handle = gzip.open(fastq_file, 'rb')
            first_1000_bytes = file_handle.read(1000)
            file_handle.close()
            
            # Rough estimation: 4 lines per read
            estimated_reads = (total_size / len(first_1000_bytes)) * (count / 4) * 1000
            return int(estimated_reads)
        except Exception as e:
            print(f"Warning: Could not count reads in {fastq_file}: {e}")
            return 0
    
    def get_summary(self):
        """Get summary statistics"""
        total_samples = len(self.samples)
        total_size = sum(s['total_size_gb'] for s in self.samples.values())
        avg_size = total_size / total_samples if total_samples > 0 else 0
        
        return {
            'total_samples': total_samples,
            'total_size_gb': total_size,
            'avg_size_gb': avg_size,
            'generation_date': datetime.now().strftime("%B %d, %Y")
        }


class ReportGenerator:
    """Generate PDF report"""
    
    def __init__(self, output_path, author, clc_version, library_kit, fastq_dir):
        self.output_path = output_path
        self.author = author
        self.clc_version = clc_version
        self.library_kit = library_kit
        self.extractor = FASTQMetadataExtractor(fastq_dir)
        self.summary = self.extractor.get_summary()
    
    def generate(self):
        """Generate the PDF report"""
        doc = SimpleDocTemplate(
            self.output_path,
            pagesize=letter,
            rightMargin=0.75*inch,
            leftMargin=0.75*inch,
            topMargin=0.75*inch,
            bottomMargin=0.75*inch,
        )
        
        elements = []
        styles = getSampleStyleSheet()
        
        # Custom styles
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            textColor=colors.HexColor('#1f4788'),
            spaceAfter=12,
            alignment=TA_CENTER,
            fontName='Helvetica-Bold'
        )
        
        heading_style = ParagraphStyle(
            'CustomHeading',
            parent=styles['Heading2'],
            fontSize=14,
            textColor=colors.HexColor('#1f4788'),
            spaceAfter=12,
            spaceBefore=12,
            fontName='Helvetica-Bold'
        )
        
        # Title
        elements.append(Paragraph("Small RNA Sequencing Project Report", title_style))
        elements.append(Spacer(1, 0.2*inch))
        
        # Project Information
        elements.append(Paragraph("Project Information", heading_style))
        project_info = [
            ['Generation Date:', self.summary['generation_date']],
            ['Author:', self.author],
            ['CLC Genomics Workbench Version:', self.clc_version],
            ['Library Kit:', self.library_kit],
        ]
        
        info_table = Table(project_info, colWidths=[2.7*inch, 3.3*inch])
        info_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (0, -1), colors.HexColor('#E8EEF7')),
            ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
        ]))
        elements.append(info_table)
        elements.append(Spacer(1, 0.3*inch))
        
        # Summary Statistics
        elements.append(Paragraph("Summary Statistics", heading_style))
        summary_data = [
            ['Metric', 'Value'],
            ['Total Samples', str(self.summary['total_samples'])],
            ['Total Data Size', f"{self.summary['total_size_gb']:.2f} GB"],
            ['Average Size per Sample', f"{self.summary['avg_size_gb']:.2f} GB"],
        ]
        
        summary_table = Table(summary_data, colWidths=[2.5*inch, 2.5*inch])
        summary_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
            ('TOPPADDING', (0, 0), (-1, -1), 8),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F0F0F0')]),
        ]))
        elements.append(summary_table)
        elements.append(Spacer(1, 0.3*inch))
        
        # General Pipeline
        elements.append(Paragraph("General Pipeline", heading_style))
        
        body_style = ParagraphStyle(
            'CustomBody',
            parent=styles['BodyText'],
            fontSize=10,
            alignment=TA_LEFT,
            spaceAfter=10,
            leading=14
        )
        
        # FastQC section
        elements.append(Paragraph("<b>FastQC v0.11.9</b>", styles['Heading3']))
        fastqc_text = (
            "FastQC was run on all FASTQ files to determine the sequencing data quality and presence of adapters. "
            "The HTML reports were exported and analyzed for quality metrics."
        )
        elements.append(Paragraph(fastqc_text, body_style))
        elements.append(Spacer(1, 0.1*inch))
        
        # Cutadapt section
        elements.append(Paragraph("<b>Cutadapt v2.10</b>", styles['Heading3']))
        cutadapt_text = (
            "Cutadapt was used to trim the reads according to Revvity's NEXTflex Small RNA v4 trimming guidelines. "
            "Due to the small size of small RNA inserts, only Read1 was used for analysis. The Nextera transposase "
            "sequence was first trimmed from Read1. Then the 3' adapter sequence (TGGAATTCTCGGGTGCCAAGG) was trimmed, "
            "while removing reads that are less than 15bp. Reads with length-based quality cutoff at 20 were retained."
        )
        elements.append(Paragraph(cutadapt_text, body_style))
        elements.append(Spacer(1, 0.1*inch))
        
        # CLC Genomics section
        elements.append(Paragraph(f"<b>CLC Genomics Workbench {self.clc_version}</b>", styles['Heading3']))
        clc_text = (
            "The trimmed reads were imported into QIAGEN CLC Genomics Workbench for further processing. "
            "CLC 'Trim Reads' was first used to remove any reads whose insert was greater than 55bp. "
            "Then 'Quantify miRNA v1.3' was used to quantify and count the miRNA using miRBase-Release_v22, "
            "prioritizing Homo sapiens with the following parameters: allow length-based isomiRs, 2 base additional "
            "up/downstream, 2 base missing up/downstream, maximum mismatch = 2, minimum sequence length = 15, "
            "maximum sequence length = 55, and minimum supporting count = 1. In addition, 'Combined miRNA Report v1.0' "
            "was used to consolidate all individual miRNA quantification reports into one report. The report was exported "
            "as an pdf file. The individual miRNA quantification reports and results were exported and zipped in excel format."
        )
        elements.append(Paragraph(clc_text, body_style))
        elements.append(Spacer(1, 0.2*inch))
        
        # Sample Details
        elements.append(Paragraph("Sample Details", heading_style))
        
        sample_data = [['Sample Name', 'i7 Index', 'i5 Index', 'Est. Reads']]
        
        for sample_name in sorted(self.extractor.samples.keys()):
            sample_info = self.extractor.samples[sample_name]
            sample_data.append([
                sample_name,
                sample_info['i7_index'],
                sample_info['i5_index'],
                f"{sample_info['read_count']:,.0f}" if sample_info['read_count'] > 0 else 'N/A',
            ])
        
        samples_table = Table(sample_data, colWidths=[2.5*inch, 1.25*inch, 1.25*inch, 1.25*inch])
        samples_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 9),
            ('FONTSIZE', (0, 1), (-1, -1), 8),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
            ('TOPPADDING', (0, 0), (-1, -1), 6),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F5F5F5')]),
        ]))
        elements.append(samples_table)
        elements.append(Spacer(1, 0.3*inch))
        
        # References
        elements.append(PageBreak())
        elements.append(Paragraph("References", heading_style))
        
        reference_style = ParagraphStyle(
            'Reference',
            parent=styles['BodyText'],
            fontSize=9,
            leftIndent=0.2*inch,
            spaceAfter=8,
            leading=11,
            textColor=colors.black
        )
        
        references = [
            "Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc",
            "Martin, Marcel. \"Cutadapt removes adapter sequences from high-throughput sequencing reads.\" EMBnet.journal [Online], 17.1 (2011): pp. 10-12. Web. 22 Jul. 2024",
            "QIAGEN CLC Main Workbench 25.0. https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/2500/index.php",
        ]
        
        for i, ref in enumerate(references, 1):
            elements.append(Paragraph(f"<b>{i}.</b> {ref}", reference_style))
        
        # Build PDF
        doc.build(elements)
        print(f"Report generated successfully: {self.output_path}")


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Generate a small RNA sequencing project report'
    )
    parser.add_argument(
        '--fastq-dir',
        default='data/FASTQ',
        help='Path to FASTQ directory (default: data/FASTQ)'
    )
    parser.add_argument(
        '--output',
        default='SmallRNA_Project_Report.pdf',
        help='Output PDF path (default: SmallRNA_Project_Report.pdf)'
    )
    parser.add_argument(
        '--author',
        default='Kevin Stachelek',
        help='Report author name'
    )
    parser.add_argument(
        '--clc-version',
        default='v25',
        help='CLC Genomics Workbench version'
    )
    parser.add_argument(
        '--library-kit',
        default='Revvity NEXTflex Small RNA v4',
        help='Library kit used for sequencing'
    )
    
    args = parser.parse_args()
    
    # Validate FASTQ directory
    if not os.path.isdir(args.fastq_dir):
        print(f"Error: FASTQ directory not found: {args.fastq_dir}")
        return 1
    
    # Generate report
    generator = ReportGenerator(
        output_path=args.output,
        author=args.author,
        clc_version=args.clc_version,
        library_kit=args.library_kit,
        fastq_dir=args.fastq_dir
    )
    
    generator.generate()
    return 0


if __name__ == '__main__':
    exit(main())
