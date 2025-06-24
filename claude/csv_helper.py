#!/usr/bin/env python3
"""
CSV Helper Script for Memory-Efficient CSV Operations
Designed to minimize token usage when working with CSV files in Claude sessions
"""

import argparse
import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
from typing import List, Optional, Dict, Any
import json


class CSVHelper:
    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        if not self.filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        # Get file info
        self.file_size = self.filepath.stat().st_size
        self.file_size_mb = self.file_size / (1024 * 1024)
        
        # Get line count efficiently
        with open(self.filepath, 'r') as f:
            self.line_count = sum(1 for _ in f)
    
    def peek(self, rows: int = 10) -> None:
        """Peek at file structure without loading entire file"""
        print(f"File: {self.filepath}")
        print(f"Size: {self.file_size_mb:.2f} MB")
        print(f"Lines: {self.line_count:,}")
        print("-" * 50)
        
        # Read sample
        df = pd.read_csv(self.filepath, nrows=rows)
        
        print(f"Columns ({len(df.columns)}): {', '.join(df.columns)}")
        print(f"Shape: {df.shape}")
        print("\nData types:")
        print(df.dtypes)
        print("\nFirst {rows} rows:")
        print(df)
    
    def read(self, start: Optional[int] = None, count: Optional[int] = None, 
             columns: Optional[List[str]] = None) -> pd.DataFrame:
        """Read specific portions of the CSV"""
        kwargs = {}
        
        if columns:
            kwargs['usecols'] = columns
        
        if start is not None:
            if count is not None:
                kwargs['skiprows'] = range(1, start)
                kwargs['nrows'] = count
            else:
                kwargs['skiprows'] = range(1, start)
        elif count is not None:
            kwargs['nrows'] = count
        
        df = pd.read_csv(self.filepath, **kwargs)
        print(f"Read {len(df)} rows, {len(df.columns)} columns")
        return df
    
    def stats(self, columns: Optional[List[str]] = None, sample_size: int = 10000) -> Dict[str, Any]:
        """Get statistics without loading entire file"""
        # Use sampling for large files
        if self.line_count > sample_size * 2:
            skip_rows = sorted(np.random.choice(range(1, self.line_count), 
                                              self.line_count - sample_size - 1, 
                                              replace=False))
            df = pd.read_csv(self.filepath, skiprows=skip_rows, usecols=columns)
            print(f"Using sample of {len(df)} rows")
        else:
            df = pd.read_csv(self.filepath, usecols=columns)
        
        stats = {
            'shape': df.shape,
            'columns': df.columns.tolist(),
            'dtypes': df.dtypes.to_dict(),
            'missing': df.isnull().sum().to_dict(),
            'numeric_stats': {}
        }
        
        # Get numeric column statistics
        numeric_cols = df.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            stats['numeric_stats'] = df[numeric_cols].describe().to_dict()
        
        return stats
    
    def missing(self, sample_size: int = 10000) -> pd.DataFrame:
        """Find missing values efficiently"""
        # Sample for large files
        if self.line_count > sample_size * 2:
            df = pd.read_csv(self.filepath, nrows=sample_size)
            print(f"Analyzing sample of {sample_size} rows")
        else:
            df = pd.read_csv(self.filepath)
        
        missing_counts = df.isnull().sum()
        missing_pct = (missing_counts / len(df)) * 100
        
        result = pd.DataFrame({
            'column': missing_counts.index,
            'missing_count': missing_counts.values,
            'missing_pct': missing_pct.values
        })
        
        return result[result['missing_count'] > 0].sort_values('missing_count', ascending=False)
    
    def duplicates(self, columns: Optional[List[str]] = None) -> int:
        """Count duplicates efficiently using chunks"""
        seen = set()
        duplicate_count = 0
        
        for chunk in pd.read_csv(self.filepath, chunksize=10000, usecols=columns):
            if columns:
                # Check specific columns
                for _, row in chunk.iterrows():
                    key = tuple(row[columns].values)
                    if key in seen:
                        duplicate_count += 1
                    else:
                        seen.add(key)
            else:
                # Check full rows
                chunk_dupes = chunk.duplicated().sum()
                duplicate_count += chunk_dupes
        
        return duplicate_count
    
    def outliers(self, column: str, method: str = 'iqr') -> pd.DataFrame:
        """Find outliers in a numeric column"""
        # Read only the needed column
        df = pd.read_csv(self.filepath, usecols=[column])
        
        if method == 'iqr':
            Q1 = df[column].quantile(0.25)
            Q3 = df[column].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            outliers = df[(df[column] < lower_bound) | (df[column] > upper_bound)]
        elif method == 'zscore':
            z_scores = np.abs((df[column] - df[column].mean()) / df[column].std())
            outliers = df[z_scores > 3]
        else:
            raise ValueError(f"Unknown method: {method}")
        
        return outliers
    
    def filter(self, column: str, value: Any, output_file: Optional[str] = None) -> pd.DataFrame:
        """Filter rows matching criteria"""
        results = []
        
        for chunk in pd.read_csv(self.filepath, chunksize=10000):
            filtered = chunk[chunk[column] == value]
            results.append(filtered)
        
        result_df = pd.concat(results, ignore_index=True)
        
        if output_file:
            result_df.to_csv(output_file, index=False)
            print(f"Saved {len(result_df)} rows to {output_file}")
        
        return result_df
    
    def fix(self, output_file: str, fix_missing: bool = False, 
            fix_types: bool = False, remove_duplicates: bool = False) -> None:
        """Fix common issues in chunks"""
        first_chunk = True
        
        for chunk in pd.read_csv(self.filepath, chunksize=10000):
            if fix_missing:
                # Fix missing values with appropriate defaults
                for col in chunk.columns:
                    if chunk[col].dtype == 'object':
                        chunk[col].fillna('', inplace=True)
                    else:
                        chunk[col].fillna(0, inplace=True)
            
            if fix_types:
                # Try to infer better types
                chunk = chunk.infer_objects()
            
            if remove_duplicates:
                chunk = chunk.drop_duplicates()
            
            # Write chunk
            if first_chunk:
                chunk.to_csv(output_file, index=False, mode='w')
                first_chunk = False
            else:
                chunk.to_csv(output_file, index=False, mode='a', header=False)
        
        print(f"Fixed file saved to {output_file}")
    
    def schema(self) -> Dict[str, Any]:
        """Get comprehensive table schema and context"""
        # Read first chunk to understand structure
        df_sample = pd.read_csv(self.filepath, nrows=100)
        
        # Try to detect index column
        first_col = df_sample.columns[0]
        has_index = df_sample[first_col].is_unique
        
        # Check if it's not just row numbers
        if has_index and df_sample[first_col].dtype == 'object':
            try:
                has_index = not df_sample[first_col].str.match(r'^\d+$').all()
            except:
                pass
        
        schema = {
            'file': str(self.filepath),
            'file_size_mb': self.file_size_mb,
            'total_rows': self.line_count - 1,  # Excluding header
            'total_columns': len(df_sample.columns),
            'columns': df_sample.columns.tolist(),
            'dtypes': {col: str(dtype) for col, dtype in df_sample.dtypes.items()},
            'index_column': first_col if has_index else None,
            'memory_usage_mb': df_sample.memory_usage(deep=True).sum() / 1024 / 1024 * (self.line_count / 100),
            'sample_values': {}
        }
        
        # Add sample values for each column
        for col in df_sample.columns[:10]:  # Limit to first 10 columns
            unique_values = df_sample[col].dropna().unique()
            if len(unique_values) > 5:
                schema['sample_values'][col] = unique_values[:5].tolist() + ['...']
            else:
                schema['sample_values'][col] = unique_values.tolist()
        
        return schema
    
    def indices(self, sample_size: int = 20) -> Dict[str, Any]:
        """Analyze table indices"""
        df_sample = pd.read_csv(self.filepath, nrows=min(sample_size, self.line_count - 1))
        
        # Check if first column is likely an index
        first_col = df_sample.columns[0]
        
        result = {
            'index_column': first_col,
            'index_type': str(df_sample[first_col].dtype),
            'is_unique': df_sample[first_col].is_unique,
            'sample_indices': df_sample[first_col].tolist(),
            'total_indices': self.line_count - 1
        }
        
        # Check if numeric or string index
        if pd.api.types.is_numeric_dtype(df_sample[first_col]):
            result['numeric_range'] = {
                'min': df_sample[first_col].min(),
                'max': df_sample[first_col].max()
            }
        
        return result
    
    def get_row(self, index_value: Optional[str] = None, position: Optional[int] = None, 
                index_col: int = 0) -> Dict[str, Any]:
        """Get single row as dictionary by index value or position"""
        if index_value is not None:
            # Read with index column
            df = pd.read_csv(self.filepath, index_col=index_col)
            if index_value in df.index:
                return df.loc[index_value].to_dict()
            else:
                raise ValueError(f"Index '{index_value}' not found")
        elif position is not None:
            # Read specific row by position
            df = pd.read_csv(self.filepath, skiprows=range(1, position + 1), nrows=1)
            if len(df) > 0:
                return df.iloc[0].to_dict()
            else:
                raise ValueError(f"Position {position} out of range")
        else:
            raise ValueError("Must specify either index_value or position")
    
    def get_rows(self, indices: List[str], index_col: int = 0) -> pd.DataFrame:
        """Get multiple rows by index values"""
        df = pd.read_csv(self.filepath, index_col=index_col)
        found_indices = [idx for idx in indices if idx in df.index]
        if not found_indices:
            raise ValueError(f"None of the indices {indices} were found")
        return df.loc[found_indices]
    
    def get_cell(self, index_value: str, column: str, index_col: int = 0) -> Any:
        """Get specific cell value by index and column"""
        df = pd.read_csv(self.filepath, index_col=index_col)
        if index_value not in df.index:
            raise ValueError(f"Index '{index_value}' not found")
        if column not in df.columns:
            raise ValueError(f"Column '{column}' not found")
        return df.loc[index_value, column]
    
    def get_cells(self, indices: List[str], columns: List[str], index_col: int = 0) -> pd.DataFrame:
        """Get multiple cells by indices and columns"""
        df = pd.read_csv(self.filepath, index_col=index_col)
        found_indices = [idx for idx in indices if idx in df.index]
        found_columns = [col for col in columns if col in df.columns]
        if not found_indices or not found_columns:
            raise ValueError("Some indices or columns not found")
        return df.loc[found_indices, found_columns]
    
    def search_index(self, pattern: Optional[str] = None, 
                    min_val: Optional[float] = None, 
                    max_val: Optional[float] = None,
                    index_col: int = 0) -> List[str]:
        """Search for indices matching pattern or range"""
        # Read only index column
        df = pd.read_csv(self.filepath, usecols=[index_col])
        index_series = df.iloc[:, 0]
        
        if pattern:
            # String pattern matching
            matches = index_series[index_series.astype(str).str.contains(pattern, regex=True)]
            return matches.tolist()
        elif min_val is not None or max_val is not None:
            # Numeric range matching
            numeric_indices = pd.to_numeric(index_series, errors='coerce')
            mask = pd.Series([True] * len(numeric_indices))
            if min_val is not None:
                mask &= numeric_indices >= min_val
            if max_val is not None:
                mask &= numeric_indices <= max_val
            return index_series[mask].tolist()
        else:
            raise ValueError("Must specify either pattern or numeric range")
    
    def filter_row(self, index_value: str, max_value: Optional[float] = None, 
                  min_value: Optional[float] = None, index_col: int = 0) -> Dict[str, Any]:
        """Filter a single row's values by numeric thresholds"""
        row_dict = self.get_row(index_value=index_value, index_col=index_col)
        
        filtered = {}
        for col, val in row_dict.items():
            try:
                numeric_val = float(val)
                if max_value is not None and numeric_val > max_value:
                    continue
                if min_value is not None and numeric_val < min_value:
                    continue
                filtered[col] = val
            except (ValueError, TypeError):
                # Keep non-numeric values
                filtered[col] = val
        
        return filtered
    
    def column_frequencies(self, columns: List[str], top_n: int = 10, 
                          exclude_values: Optional[List[str]] = None,
                          residue_format: bool = False) -> Dict[str, Dict[str, Any]]:
        """Get value frequencies for specified columns
        
        Args:
            columns: List of column names to analyze
            top_n: Number of top values to return
            exclude_values: Values to exclude from analysis
            residue_format: If True, parse values as residue+position (e.g., K296 -> K)
        """
        if exclude_values is None:
            exclude_values = ['-', 'nan', 'NaN', 'None', '']
        
        # First get the index column name
        df_header = pd.read_csv(self.filepath, nrows=0)
        index_col_name = df_header.columns[0]
        
        # Read with the specified columns plus index
        df = pd.read_csv(self.filepath, usecols=[index_col_name] + columns, index_col=0)
        total_rows = len(df)
        
        result = {}
        for col in columns:
            if col in df.columns:
                if residue_format:
                    # Extract just the residue letter from format like 'K296'
                    residue_series = df[col].apply(lambda x: x[0] if pd.notna(x) and x not in exclude_values and len(x) > 0 else None)
                    # Get value counts on residue letters only
                    counts = residue_series.value_counts()
                    # Remove None values
                    counts = counts[counts.index.notna()]
                else:
                    # Get value counts on full values
                    counts = df[col].value_counts()
                    # Filter out excluded values
                    counts = counts[~counts.index.isin(exclude_values)]
                
                # Calculate percentages
                percentages = (counts / total_rows * 100).round(1)
                
                # Get top N values
                top_values = []
                for value, count in counts.head(top_n).items():
                    top_values.append({
                        'value': str(value),
                        'count': int(count),
                        'percentage': float(percentages[value])
                    })
                
                # For residue format, also show example positions
                if residue_format and len(top_values) > 0:
                    # Get examples of positions for each residue type
                    for i, item in enumerate(top_values):
                        residue = item['value']
                        # Find all values starting with this residue
                        mask = df[col].str.startswith(residue, na=False)
                        examples = df.loc[mask, col].value_counts().head(3).index.tolist()
                        item['position_examples'] = examples
                
                result[col] = {
                    'total_rows': total_rows,
                    'unique_values': len(counts),
                    'top_values': top_values
                }
        
        return result
    
    def report(self, output_file: Optional[str] = None, sample_size: int = 10000) -> str:
        """Generate comprehensive report"""
        report_lines = []
        report_lines.append(f"CSV Analysis Report: {self.filepath}")
        report_lines.append("=" * 60)
        report_lines.append(f"File Size: {self.file_size_mb:.2f} MB")
        report_lines.append(f"Total Rows: {self.line_count:,}")
        
        # Get statistics
        stats = self.stats(sample_size=sample_size)
        report_lines.append(f"Total Columns: {len(stats['columns'])}")
        report_lines.append("\nColumn Information:")
        report_lines.append("-" * 40)
        
        for col, dtype in stats['dtypes'].items():
            missing = stats['missing'].get(col, 0)
            report_lines.append(f"{col}: {dtype} (missing: {missing})")
        
        # Missing value summary
        total_missing = sum(stats['missing'].values())
        if total_missing > 0:
            report_lines.append(f"\nTotal Missing Values: {total_missing}")
            report_lines.append("Columns with missing values:")
            for col, count in sorted(stats['missing'].items(), key=lambda x: x[1], reverse=True):
                if count > 0:
                    pct = (count / stats['shape'][0]) * 100
                    report_lines.append(f"  {col}: {count} ({pct:.1f}%)")
        
        # Numeric statistics
        if stats['numeric_stats']:
            report_lines.append("\nNumeric Column Statistics:")
            report_lines.append("-" * 40)
            for col, col_stats in stats['numeric_stats'].items():
                report_lines.append(f"\n{col}:")
                for stat, value in col_stats.items():
                    report_lines.append(f"  {stat}: {value:.2f}")
        
        report = "\n".join(report_lines)
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(report)
            print(f"Report saved to {output_file}")
        else:
            print(report)
        
        return report


def main():
    parser = argparse.ArgumentParser(description='Memory-efficient CSV operations')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Peek command
    peek_parser = subparsers.add_parser('peek', help='Peek at file structure')
    peek_parser.add_argument('file', help='CSV file path')
    peek_parser.add_argument('--rows', type=int, default=10, help='Number of rows to show')
    
    # Read command
    read_parser = subparsers.add_parser('read', help='Read specific portions')
    read_parser.add_argument('file', help='CSV file path')
    read_parser.add_argument('--start', type=int, help='Starting row (1-indexed)')
    read_parser.add_argument('--count', type=int, help='Number of rows to read')
    read_parser.add_argument('--columns', help='Comma-separated column names')
    
    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Get statistics')
    stats_parser.add_argument('file', help='CSV file path')
    stats_parser.add_argument('--columns', help='Comma-separated column names')
    stats_parser.add_argument('--sample', type=int, default=10000, help='Sample size')
    
    # Missing command
    missing_parser = subparsers.add_parser('missing', help='Find missing values')
    missing_parser.add_argument('file', help='CSV file path')
    missing_parser.add_argument('--sample', type=int, default=10000, help='Sample size')
    
    # Duplicates command
    dup_parser = subparsers.add_parser('duplicates', help='Count duplicates')
    dup_parser.add_argument('file', help='CSV file path')
    dup_parser.add_argument('--columns', help='Comma-separated column names')
    
    # Outliers command
    outlier_parser = subparsers.add_parser('outliers', help='Find outliers')
    outlier_parser.add_argument('file', help='CSV file path')
    outlier_parser.add_argument('--column', required=True, help='Column to analyze')
    outlier_parser.add_argument('--method', choices=['iqr', 'zscore'], default='iqr')
    
    # Filter command
    filter_parser = subparsers.add_parser('filter', help='Filter rows')
    filter_parser.add_argument('file', help='CSV file path')
    filter_parser.add_argument('--column', required=True, help='Column to filter on')
    filter_parser.add_argument('--value', required=True, help='Value to match')
    filter_parser.add_argument('--output', help='Output file path')
    
    # Fix command
    fix_parser = subparsers.add_parser('fix', help='Fix common issues')
    fix_parser.add_argument('file', help='Input CSV file path')
    fix_parser.add_argument('output', help='Output CSV file path')
    fix_parser.add_argument('--fix-missing', action='store_true', help='Fill missing values')
    fix_parser.add_argument('--fix-types', action='store_true', help='Fix data types')
    fix_parser.add_argument('--remove-duplicates', action='store_true', help='Remove duplicates')
    
    # Report command
    report_parser = subparsers.add_parser('report', help='Generate report')
    report_parser.add_argument('file', help='CSV file path')
    report_parser.add_argument('--output', help='Output file path')
    report_parser.add_argument('--sample', type=int, default=10000, help='Sample size')
    
    # Schema command
    schema_parser = subparsers.add_parser('schema', help='Get table schema and context')
    schema_parser.add_argument('file', help='CSV file path')
    
    # Indices command
    indices_parser = subparsers.add_parser('indices', help='Analyze table indices')
    indices_parser.add_argument('file', help='CSV file path')
    indices_parser.add_argument('--sample', type=int, default=20, help='Number of indices to sample')
    
    # Get-row command
    getrow_parser = subparsers.add_parser('get-row', help='Get single row as dictionary')
    getrow_parser.add_argument('file', help='CSV file path')
    getrow_parser.add_argument('--index', help='Index value to retrieve')
    getrow_parser.add_argument('--position', type=int, help='Row position (0-based)')
    getrow_parser.add_argument('--index-col', type=int, default=0, help='Index column position')
    
    # Get-rows command
    getrows_parser = subparsers.add_parser('get-rows', help='Get multiple rows')
    getrows_parser.add_argument('file', help='CSV file path')
    getrows_parser.add_argument('--indices', required=True, help='Comma-separated index values')
    getrows_parser.add_argument('--index-col', type=int, default=0, help='Index column position')
    
    # Get-cell command
    getcell_parser = subparsers.add_parser('get-cell', help='Get specific cell value')
    getcell_parser.add_argument('file', help='CSV file path')
    getcell_parser.add_argument('--index', required=True, help='Row index value')
    getcell_parser.add_argument('--column', required=True, help='Column name')
    getcell_parser.add_argument('--index-col', type=int, default=0, help='Index column position')
    
    # Get-cells command
    getcells_parser = subparsers.add_parser('get-cells', help='Get multiple cells')
    getcells_parser.add_argument('file', help='CSV file path')
    getcells_parser.add_argument('--indices', required=True, help='Comma-separated index values')
    getcells_parser.add_argument('--columns', required=True, help='Comma-separated column names')
    getcells_parser.add_argument('--index-col', type=int, default=0, help='Index column position')
    
    # Search-index command
    searchindex_parser = subparsers.add_parser('search-index', help='Search indices')
    searchindex_parser.add_argument('file', help='CSV file path')
    searchindex_parser.add_argument('--pattern', help='Regex pattern to match')
    searchindex_parser.add_argument('--min', type=float, help='Minimum value for numeric indices')
    searchindex_parser.add_argument('--max', type=float, help='Maximum value for numeric indices')
    searchindex_parser.add_argument('--index-col', type=int, default=0, help='Index column position')
    
    # Filter-row command
    filterrow_parser = subparsers.add_parser('filter-row', help='Filter row values by threshold')
    filterrow_parser.add_argument('file', help='CSV file path')
    filterrow_parser.add_argument('--index', required=True, help='Row index value')
    filterrow_parser.add_argument('--min-value', type=float, help='Minimum value threshold')
    filterrow_parser.add_argument('--max-value', type=float, help='Maximum value threshold')
    filterrow_parser.add_argument('--index-col', type=int, default=0, help='Index column position')
    
    # Column-frequencies command
    colfreq_parser = subparsers.add_parser('column-frequencies', help='Get value frequencies for columns')
    colfreq_parser.add_argument('file', help='CSV file path')
    colfreq_parser.add_argument('--columns', required=True, help='Comma-separated column names')
    colfreq_parser.add_argument('--top-n', type=int, default=10, help='Number of top values to show')
    colfreq_parser.add_argument('--exclude', help='Comma-separated values to exclude (default: -,nan,NaN,None,empty)')
    colfreq_parser.add_argument('--residue-format', action='store_true', help='Parse values as residue+position format (e.g., K296 -> K)')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    try:
        helper = CSVHelper(args.file)
        
        if args.command == 'peek':
            helper.peek(args.rows)
        
        elif args.command == 'read':
            columns = args.columns.split(',') if args.columns else None
            df = helper.read(args.start, args.count, columns)
            print(df)
        
        elif args.command == 'stats':
            columns = args.columns.split(',') if args.columns else None
            stats = helper.stats(columns, args.sample)
            print(json.dumps(stats, indent=2, default=str))
        
        elif args.command == 'missing':
            result = helper.missing(args.sample)
            print(result)
        
        elif args.command == 'duplicates':
            columns = args.columns.split(',') if args.columns else None
            count = helper.duplicates(columns)
            print(f"Found {count} duplicate rows")
        
        elif args.command == 'outliers':
            outliers = helper.outliers(args.column, args.method)
            print(f"Found {len(outliers)} outliers:")
            print(outliers)
        
        elif args.command == 'filter':
            result = helper.filter(args.column, args.value, args.output)
            if not args.output:
                print(result)
        
        elif args.command == 'fix':
            helper.fix(args.output, args.fix_missing, args.fix_types, args.remove_duplicates)
        
        elif args.command == 'report':
            helper.report(args.output, args.sample)
        
        elif args.command == 'schema':
            schema = helper.schema()
            print(json.dumps(schema, indent=2, default=str))
        
        elif args.command == 'indices':
            indices = helper.indices(args.sample)
            print(json.dumps(indices, indent=2, default=str))
        
        elif args.command == 'get-row':
            row = helper.get_row(args.index, args.position, args.index_col)
            print(json.dumps(row, indent=2, default=str))
        
        elif args.command == 'get-rows':
            indices = args.indices.split(',')
            df = helper.get_rows(indices, args.index_col)
            print(df.to_json(orient='index', indent=2))
        
        elif args.command == 'get-cell':
            value = helper.get_cell(args.index, args.column, args.index_col)
            print(f"{value}")
        
        elif args.command == 'get-cells':
            indices = args.indices.split(',')
            columns = args.columns.split(',')
            df = helper.get_cells(indices, columns, args.index_col)
            print(df.to_json(orient='index', indent=2))
        
        elif args.command == 'search-index':
            results = helper.search_index(args.pattern, args.min, args.max, args.index_col)
            print(f"Found {len(results)} matching indices:")
            for idx in results:
                print(f"  {idx}")
        
        elif args.command == 'filter-row':
            filtered = helper.filter_row(args.index, args.max_value, args.min_value, args.index_col)
            print(json.dumps(filtered, indent=2, default=str))
        
        elif args.command == 'column-frequencies':
            columns = args.columns.split(',')
            exclude_values = args.exclude.split(',') if args.exclude else None
            result = helper.column_frequencies(columns, args.top_n, exclude_values, args.residue_format)
            print(json.dumps(result, indent=2, default=str))
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()