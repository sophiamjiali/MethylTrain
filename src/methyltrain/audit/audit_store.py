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
from typing import List
import pandas as pd
from datetime import datetime

class AuditStore:
    
    # ~~~~~| I/O of Schema |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def __init__(self, db_path):
        self.conn = sqlite3.connect(db_path)
        self._init_schema()

    def _init_schema(self):
        self.conn.execute(
            """
            CREATE TABLE IF NOT EXISTS audit (
                file_id TEXT PRIMARY KEY,
                file_name TEXT,

                download_status INTEGER DEFAULT 0,
                metadata_status INTEGER DEFAULT 0,
                biospecimen_status INTEGER DEFAULT 0,
                qc_status INTEGER DEFAULT 0,

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
            """
            INSERT OR IGNORE INTO audit (file_id, updated_at)
            VALUES (?, ?)
            """, [
                (fid, self._now()) for fid in file_ids
            ]
        )
        self.conn.commit()

    def to_dataframe(self) -> pd.DataFrame:
        return pd.read_sql_query("SELECT * FROM audit", self.conn)

    def export_csv(self, path: str):
        df = self.to_dataframe()
        df.to_csv(path, index=False)

    def _now(self) -> str:
        return datetime.utcnow().isoformat()
    
    def close(self):
        self.conn.close()


    # ~~~~~| Pipeline Status Updates |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    def set_download_status(self, file_id: str, file_name: str, status: int):
        self.conn.execute(
            """
            UPDATE audit
            SET file_name = ?,
                download_status = ?,
                download_timestamp = ?,
                updated_at = ?
            WHERE file_id = ?
            """, (file_name, status, self._now(), self._now(), file_id)
        )
        self.conn.commit()

    def set_metadata_status(self, file_id: str, status: int):
        self.conn.execute(
            """
            UPDATE audit
            SET metadata_status = ?,
                metadata_timestamp = ?,
                updated_at = ?
            WHERE file_id = ?
            """, (status, self._now(), self._now(), file_id)
        )
        self.conn.commit()

    def set_biospecimen_status(self, file_id: str, status: int):
        self.conn.execute(
            """
            UPDATE audit
            SET biospecimen_status = ?,
                biospecimen_timestamp = ?,
                updated_at = ?
            WHERE file_id = ?
            """, (status, self._now(), self._now(), file_id)
        )
        self.conn.commit()

    def set_qc_status(self, file_id: str, status: int):
        self.conn.execute(
            """
            UPDATE audit
            SET qc_status = ?,
                qc_timestamp = ?,
                updated_at = ?
            WHERE file_id = ?
            """, (status, self._now(), self._now(), file_id)
        )
        self.conn.commit()

    def set_raw_path(self, file_id: str, path: str):
        self.conn.execute(
            """
            UPDATE audit
            SET raw_data_path = ?,
                updated_at = ?
            WHERE file_id = ?
            """, (path, self._now(), file_id)
        )
        self.conn.commit()

    # ~~~~~| Getters |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def get_ids_by_download_status(self, status: int) -> list:
        cursor = self.conn.execute(
            """
            SELECT file_id
            FROM audit
            WHERE download_status = ?
            """, (status,)
        )
        return [row[0] for row in cursor.fetchall()]

# [END]