# ==============================================================================
# Script:           audit_store.py
# Purpose:          Handles audit store initializing and logging
# Author:           Sophia Li
# Affiliation:      CCG Lab, Princess Margaret Cancer Center, UHN, UofT
# Date:             2026-01-21
#
# Notes:            Maintains the audit store as an SQL file for live updating,
#                   exporting it as a CSV upon completion.
# ==============================================================================

import sqlite3

import pandas as pd

from typing import List
from datetime import datetime
from pathlib import Path

from methyltrain.constants.database import AUDIT_FIELDS, AUDIT_NAME

class AuditStore:
    
    # ~~~~~| I/O of Schema |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, db_path):
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self._init_schema()

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc, tb):
        self.close()

    def _init_schema(self):
        self.conn.execute(
            f"""
            CREATE TABLE IF NOT EXISTS {AUDIT_NAME} (
                file_id TEXT PRIMARY KEY,
                file_name TEXT,

                download_status INTEGER DEFAULT 0,
                metadata_status INTEGER DEFAULT 0,
                biospecimen_status INTEGER DEFAULT 0,
                qc_pass INTEGER DEFAULT 0,

                download_timestamp TEXT,
                metadata_timestamp TEXT,
                biospecimen_timestamp TEXT,
                qc_timestamp TEXT,
                
                raw_data_path TEXT,
                updated_at TEXT
            );
            """
        )
        self.conn.commit()


    def initialize(self, file_ids: List[str]):
        self.conn.executemany(
            f"""
            INSERT OR IGNORE INTO {AUDIT_NAME} (file_id, updated_at)
            VALUES (?, ?)
            """, [
                (fid, self._now()) for fid in file_ids
            ]
        )
        self.conn.commit()

    def to_dataframe(self) -> pd.DataFrame:
        return pd.read_sql_query(f"SELECT * FROM {AUDIT_NAME}", self.conn)

    def export_csv(self, path: Path):
        df = self.to_dataframe()
        df.to_csv(path, index=False)

    def _now(self) -> str:
        return datetime.utcnow().isoformat()
    
    def close(self):
        self.conn.close()


    # ~~~~~| Pipeline Status Updates |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def apply_download_report(self, report: list[dict]):
        """
        Updates the AuditStore with the DNA Methylation data download report.
        """
        self._apply_report('download', report)
        rows = [{
            'file_id': r['file_id'],
            'file_name': r.get('file_name')
        } for r in report]

        with self.conn:
            self.conn.executemany(
                f"""
                UPDATE {AUDIT_NAME}
                SET file_name = :file_name
                WHERE file_id = :file_id
                """, rows
            )
        return
    
    def apply_metadata_report(self, report: list[dict]):
        """
        Updates the AuditStore with the metadata download report.
        """
        self._apply_report('metadata', report)
        return
    

    def apply_biospecimen_report(self, report: list[dict]):
        """
        Updates the AuditStore with the metadata download report.
        """
        self._apply_report('biospecimen', report)
        return
    

    def apply_cleaning_report(self, report: list[dict]):
        """
        Updates the AuditStore with the cleaning download report.
        """
        now = self._now()
        rows = [{
            'file_id': r['file_id'],
            'raw_data_path': r['raw_data_path'],
            'updated_at': now,
        } for r in report]

        with self.conn:
            self.conn.executemany(
                f"""
                UPDATE {AUDIT_NAME}
                SET raw_data_path = :raw_data_path,
                    updated_at = :updated_at
                WHERE file_id = :file_id
                """, rows
            )
        return
    
    def apply_qc_report(self, report: list[dict]):
        """
        Updates the AuditStore with the quality control report.
        """
        now = self._now()
        rows = [{
            'file_id': r['file_id'],
            'qc_pass': r['qc_pass'],
            'timestamp': now,
            'updated_at': now,
        } for r in report]

        with self.conn:
            self.conn.executemany(
                f"""
                UPDATE {AUDIT_NAME}
                SET qc_pass = :qc_pass,
                    qc_timestamp = :timestamp,
                    updated_at = :updated_at
                WHERE file_id = :file_id
                """, rows
            )
        return
        

    def _apply_report(self, report_type: str, report: list[dict]):
        """
        Generic Helper for updating the AuditStore with a report.
        """
        now = self._now()
        rows = [{
            'file_id': r['file_id'],
            'status': r[f'{report_type}_status'],
            'timestamp': now,
            'updated_at': now,
        } for r in report]

        with self.conn:
            self.conn.executemany(
                f"""
                UPDATE {AUDIT_NAME}
                SET {report_type}_status = :status,
                    {report_type}_timestamp = :timestamp,
                    updated_at = :updated_at
                WHERE file_id = :file_id
                """, rows
            )
        return


    # ~~~~~| Getters |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def get_files_by_download_status(self,
                                     status: int, 
                                     columns = ("file_id",)
                                    ) -> list[tuple]:
        if not columns: 
            raise ValueError ("Atleast one column must be requested")
        
        invalid = set(columns) - AUDIT_FIELDS
        if invalid:
            raise ValueError(f"Invalid columns requested: {sorted(invalid)}")
        
        query_columns = ", ".join(columns)
        cursor = self.conn.execute(
            f"""
            SELECT {query_columns}
            FROM {AUDIT_NAME}
            WHERE download_status = ?
            """, (status,)
        )
        return cursor.fetchall()

    def get_ids_by_download_status(self, status: int) -> list:
        cursor = self.conn.execute(
            f"""
            SELECT file_id
            FROM {AUDIT_NAME}
            WHERE download_status = ?
            """, (status,)
        )
        return [row[0] for row in cursor.fetchall()]
    
# [END]